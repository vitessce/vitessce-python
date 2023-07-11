import os
from os.path import join
import requests
from jsonschema import validate, ValidationError
import scanpy as sc
import pandas as pd
import anndata
import gzip
import io
from vitessce.data_utils import (
    optimize_adata,
    VAR_CHUNK_SIZE
)


class CellBrowserToVitessceConfigConverter:

    def __init__(self, project_name, output_dir, keep_only_marker_genes):
        self.url_prefix = f"https://cells.ucsc.edu/{self.project_name}"
        self.project_name = project_name
        self.output_dir = output_dir
        self.keep_only_marker_genes = keep_only_marker_genes
        self.cellbrowser_config = {}
        self.adata = None

    def _filter_data(self):
        # Filter data
        self.adata.var_names_make_unique()
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        # Normalize the data
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

        # If marker genes are defined, keep only marker genes in the Anndata object
        if "topMarkers" in self.cellbrowser_config and len(self.cellbrowser_config["topMarkers"].keys()) > 0 and self.keep_only_marker_genes:
            print("Filtering out all non-marker genes from Anndata object ...")
            marker_genes = [gene for sublist in self.cellbrowser_config["topMarkers"].values() for gene in sublist]
            self.adata = self.adata[:, self.adata.var_names.isin(marker_genes)]
            print("Successfully filtered out all non-marker genes from Anndata object.")

    def _write_to_zarr_store(self):
        data_dir = join(self.output_dir, self.project_name)
        zarr_filepath = join(data_dir, "out.adata.zarr")

        obs_columns_list = self.adata.obs.columns.tolist()
        obsm_keys = list(self.adata.obsm.keys())
        var_cols_list = list(self.adata.var.keys())
        print("About to Anndata object to the Zarr store. The following properties will be saved:")
        print(f"  Obs columns: {obs_columns_list}")
        print(f"  Obsm keys: {obsm_keys}")
        print(f"  Var columns: {var_cols_list}")

        # This is done because of the optimize_adata function, which expects the obsm key to be numpy type
        for key in self.adata.obsm.keys():
            if isinstance(self.adata.obsm[key], pd.DataFrame):
                print(f"obsm {key} is an instance of DataFrame, converting it to numpy array.")
                self.adata.obsm[key] = self.adata.obsm[key].to_numpy()

        self.adata = optimize_adata(
            self.adata,
            obs_cols=obs_columns_list,
            obsm_keys=obsm_keys,
            optimize_X=True,
            var_cols=var_cols_list,
        )

        os.makedirs(os.path.dirname(data_dir), exist_ok=True)
        self.adata.write_zarr(zarr_filepath, chunks=[self.adata.shape[0], VAR_CHUNK_SIZE])
        print("Successfully saved Anndata object to the Zarr store.")

    def _load_expr_matrix(self):
        gzip_file = None
        expr_matrix = None

        try:
            print("Downloading expression matrix ...")
            # Collect the exprMatrix and meta url suffixes
            expr_matrix_url_suffix = self.cellbrowser_config["fileVersions"]["outMatrix"]["fname"].split("/")[-1]
            full_url = "/".join([self.url_prefix, expr_matrix_url_suffix])
            response = requests.get(full_url)
            response.raise_for_status()

            # Create a BytesIO object from the content
            gzip_file = io.BytesIO(response.content)
            print(f"Successfully downloaded expression matrix {full_url}.")
        except Exception as e:
            print(f"Could not download expression matrix {full_url} because: {e}.")
            raise e

        print("Loading expression matrix into Anndata object ...")
        with gzip.open(gzip_file, 'rt') as f:
            expr_matrix = pd.read_csv(f, sep='\t', index_col=0).T  # transpose it, because of Scanpy

        # Now create anndata object
        self.adata = anndata.AnnData(X=expr_matrix)
        self.adata.var['gene'] = self.adata.var_names
        first_gene = list(self.adata.var_names)[0]
        if len(first_gene.split("|")) == 2:
            # TODO: sometimes we might want to keep them both
            print("This dataset uses the format identifier|symbol for the ad.obs gene names (e.g. “ENSG0123123.3|HOX3”). We are keeping only the symbol.")
            self.adata.var_names = [x.split("|")[1] for x in list(self.adata.var_names)]
        print("Successfully loaded expression matrix into Anndata object.")

    def _load_coordinates(self):
        coordinate_types = {
            'tsne': 'X_tsne',
            't-sne': 'X_tsne',
            'umap': 'X_umap',
            'pca': 'X_pca',
            'spatial': 'X_spatial',
            'segmentations': 'X_segmentations'
        }

        coord_urls = {}

        for obj in self.cellbrowser_config["coords"]:
            short_label = obj["shortLabel"].lower()
            for term, key in coordinate_types.items():
                if term in short_label:
                    if 'textFname' in obj:
                        coord_urls[key] = obj['textFname']
                    else:
                        print(f"Detected {term} type of coordinates, but could not find link to coordinates files. Skipping.")
                    break

        if len(coord_urls.keys()) == 0:
            print("WARN: Dataset will be displayed without an embedding. No valid coordinates were found.")
        else:
            print(f"Successful extraction of the following coordinates and URLS: {coord_urls}")

        for (coord_type, url_suffix) in coord_urls.items():
            try:
                print(f"Adding {coord_type} to Anndata object ...")
                response = requests.get("/".join([self.url_prefix, url_suffix]))
                response.raise_for_status()
                embedding_file = io.BytesIO(response.content)

                with gzip.open(embedding_file, 'rt') as f:
                    coords = pd.read_csv(f, sep='\t', index_col=0)

                # Ensure the indices are of the same type
                coords.index = coords.index.astype(str)

                extra_in_adata = set(self.adata.obs.index) - set(coords.index)
                extra_in_coords = set(coords.index) - set(self.adata.obs.index)

                # If there's an extra cell in adata
                if extra_in_adata:
                    self.adata = self.adata[~self.adata.obs.index.isin(extra_in_adata)]

                # If there's an extra cell in coords
                if extra_in_coords:
                    coords = coords.drop(index=extra_in_coords)

                # Add the coordinates to the `obsm` attribute of `adata`
                self.adata.obsm[coord_type] = coords
                print(f"{coord_type} successfully added.")
            except Exception as e:
                print(f"Could not add {coord_type} to Anndata object because: {e}.")
                raise e
            print("Done adding coordinates to the Anndata object.")

    def _load_cell_metadata(self):
        try:
            print("Adding cell metadata to Anndata object ...")
            # Load meta
            meta_url_suffix = self.cellbrowser_config["fileVersions"]["outMeta"]["fname"].split("/")[-1]
            response = requests.get("/".join([self.url_prefix, meta_url_suffix]))
            response.raise_for_status()

            # Create a BytesIO object from the content
            meta_file = io.BytesIO(response.content)

            meta = pd.read_csv(meta_file, sep='\t', index_col=0)

            self.adata.obs = meta
            # remove space from obs column names to avoid Vitessce breaking downstream
            self.adata.obs.columns = [x.replace(" ", "") for x in self.adata.obs.columns.tolist()]
            print(f"Successfully downloaded metadata {meta_url_suffix}.")
        except Exception as e:
            print(f"Could not download metadata {meta_url_suffix} because: {e}.")
            raise e

    def _validate_config(self):
        schema = {
            "type": "object",
            "properties": {
                "fileVersions": {
                    "type": "object",
                    "properties": {
                        "outMatrix": {
                            "type": "object",
                            "properties": {
                                "fname": {"type": "string"},
                            },
                            "required": ["fname"]
                        },
                    },
                    "required": ["outMatrix"]
                },
                "coords": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "shortLabel": {"type": "string"},
                            "textFname": {"type": "string"},
                        },
                        "required": ["shortLabel", "textFname"]
                    },
                    "minItems": 1
                },
            },
            "required": ["fileVersions", "coords"]
        }

        try:
            validate(instance=self.cellbrowser_config, schema=schema)
            print("CellBrowser config is valid. Proceeding further with conversion.")
            return True
        except ValidationError as e:
            print("Invalid CellBrowser config: ", e.message)
            return False

    def _download_config(self):
        config_url = "/".join([self.url_prefix, "dataset.json"])
        try:
            response = requests.get(config_url)
            response.raise_for_status()
            self.cellbrowser_config = response.json()
            print(f"Successfully fetched configuration: {config_url}.")
            return self._validate_config()
        except Exception as e:
            print(f"Could not get configuration for dataset {self.project_name} because: {e}.")
            return False


def convert(project_name, output_dir, keep_only_marker_genes=False):
    print(f"Converting CellBrowser config for project {project_name} to Vitessce format and saving it to {output_dir}")
    config_converter = CellBrowserToVitessceConfigConverter(project_name, output_dir, keep_only_marker_genes)
    config_is_valid = config_converter._download_config()
    if config_is_valid:
        config_converter._load_expr_matrix()
        config_converter._load_cell_metadata()
        config_converter._load_coordinates()
        config_converter._filter_data()
        config_converter._write_to_zarr_store()
        print(f"CellBrowser config finished conversion. See the output files in {output_dir}.")
    else:
        print("CellBrowser config could not be converted, because it is not valid.")
