import requests
from jsonschema import validate, ValidationError
import scanpy as sc
import pandas as pd
import anndata
import gzip
import io

from vitessce.data_utils import (
    optimize_adata,
)


class CellBrowserToAnndataZarrConverter:
    """
    CellBrowserToAnndataZarrConverter class is responsible for converting a Cell Browser project to an Anndata-Zarr format.
    """

    def __init__(self, project_name, keep_only_marker_genes=False):
        """
        Constructor for the CellBrowserToAnndataZarrConverter class.
        :param project_name: The name of the Cell Browser project to be converted.
        :type project_name: str
        :param output_dir: The output directory where the converted project will be saved.
        :type output_dir: str
        :param keep_only_marker_genes: A flag indicating whether to keep only marker genes in the Anndata object.
        :type keep_only_marker_genes: bool
        """
        url_suffix = "/".join(project_name.split("+"))
        self.url_prefix = f"https://cells.ucsc.edu/{url_suffix}"
        self.project_name = project_name
        self.keep_only_marker_genes = keep_only_marker_genes
        self.cellbrowser_config = {}
        self.adata = None

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
                        "required": ["shortLabel"]
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

    def download_config(self):
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

    def _load_expr_matrix(self):
        """
        Downloads the expression matrix for the project from the CellBrowser website.
        Creates an Anndata object with it.
        """
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
            print(f"Could not download expression matrix because: {e}.")
            raise e

        print("Loading expression matrix into Anndata object ...")

        with gzip.open(gzip_file, 'rt') as f:
            expr_matrix = pd.read_csv(f, sep='\t', index_col=0).T  # transpose it, because of Scanpy

        # Now create anndata object
        self.adata = anndata.AnnData(X=expr_matrix, dtype='float32')
        self.adata.var['gene'] = self.adata.var_names

        # Filter out nan values
        self.adata = self.adata[~self.adata.obs.index.isnull(), :]
        self.adata = self.adata[:, ~self.adata.var.index.isnull()]

        first_gene = list(self.adata.var_names)[0]
        if len(first_gene.split("|")) == 2:
            # TODO: sometimes we might want to keep them both
            print("This dataset uses the format identifier|symbol for the ad.obs gene names (e.g. “ENSG0123123.3|HOX3”). We are keeping only the symbol.")
            self.adata.var_names = [x.split("|")[1] for x in list(self.adata.var_names)]

    def _load_coordinates(self):
        """
        Downloads coordinate data for the project from the CellBrowser website and loads them into the Anndata object.
        The URLs for the coordinate data are read from the cellbrowser configuration.

        Note: only coordinates that contain:
            'tsne', 't-sne', 'umap', 'pca', 'segmentations' or 'spatial'
        in their name will taken into an account.
        """
        coordinate_types = {
            'tsne': 'X_tsne',
            't-sne': 'X_tsne',
            'umap': 'X_umap',
            'pca': 'X_pca',
            'spatial': 'X_spatial',
            'segmentations': 'X_segmentations'
        }

        coord_urls = {}

        assert self.adata is not None

        for obj in self.cellbrowser_config["coords"]:
            short_label = obj["shortLabel"].lower()
            for term, key in coordinate_types.items():
                if term in short_label:
                    if 'textFname' in obj:
                        coord_urls[key] = obj['textFname']
                    else:
                        # if textFname is not defined, that means that the name of file is the shortLabel
                        labels = obj["shortLabel"].split(" ")
                        filtered_labels = [label.replace("-", "Minus") for label in labels]
                        print(filtered_labels)
                        file_name = "_".join(filtered_labels)
                        coord_urls[key] = ".".join((file_name, "coords.tsv.gz"))

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
        """
        Downloads cell metadata for the project from the CellBrowser website and loads it into the Anndata object.
        The URL for the cell metadata is read from the cellbrowser configuration.
        """
        try:
            print("Adding cell metadata to Anndata object ...")
            assert self.adata is not None
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
            print(f"Could not download metadata because: {e}.")
            raise e

    def _filter_data(self):
        """
        This function applies various preprocessing steps to the data contained in the Anndata object of the class instance.

        Firstly, it ensures the gene names in the dataset are unique. It then applies filtering to remove cells with less
        than 200 genes and genes that are expressed in less than 3 cells.

        Following the filtering, it normalizes the data by scaling the total count per cell to a target sum of 10,000.
        The values are then log-transformed for downstream analysis.

        If the flag 'keep_only_marker_genes' is set to True and marker genes are specified in the cellbrowser configuration,
        the function further filters the Anndata object to only retain these marker genes.

        Note: The preprocessing steps here, especially the filtering and normalization steps, are a common part of the
        workflow in single-cell RNA sequencing data analysis.
        They are taken from this tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/tutorial_pearson_residuals.html.
        """
        assert self.adata is not None
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

    def create_anndata_object(self):
        """
        Creates an Anndata object out of the data downloaded from the CellBrowser website. In particular, it uses:
        - The expression matrix
        - The cell metadata
        - The coordinates
        The name of each of those files is read from the cellbrowser configuration.
        """
        self._load_expr_matrix()
        self._load_cell_metadata()
        self._load_coordinates()
        self._filter_data()

    def export_anndata_object(self):
        """
        Writes the contents of the Anndata object of the class instance to a Zarr store. Saves the Zarr store to the output directory
        that the class is instantiated with.
        """
        obs_columns_list = self.adata.obs.columns.tolist()
        obsm_keys = list(self.adata.obsm.keys())
        var_cols_list = list(self.adata.var.keys())
        print("About to write the Anndata object to the Zarr store. The following properties will be saved:")
        print(f"  Obs columns: {obs_columns_list}")
        print(f"  Obsm keys: {obsm_keys}")
        print(f"  Var columns: {var_cols_list}")

        # This is done because of the optimize_adata function, which expects the obsm key to be numpy type
        for key in self.adata.obsm.keys():
            if isinstance(self.adata.obsm[key], pd.DataFrame):
                print(f"obsm {key} is an instance of DataFrame, converting it to numpy array.")
                self.adata.obsm[key] = self.adata.obsm[key].to_numpy()

        for col in self.adata.obs.columns.tolist():
            if self.adata.obs[col].dtype == object:
                self.adata.obs[col] = self.adata.obs[col].astype('category')

        return optimize_adata(
            self.adata,
            obs_cols=obs_columns_list,
            obsm_keys=obsm_keys,
            optimize_X=True,
            var_cols=var_cols_list,
        ).copy()


def convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False):
    """
    Given a CellBrowser project name, download the config, convert it to an Anndata-Zarr format,
    which is digestable by Vitessce, and save it to the output directory.

    :type project_name: str
    :param project_name: CellBrowser project name
    :param output_dir: Output directory to save the Anndata-Zarr file
    :type output_dir: str, "vitessce-files" by default
    :param keep_only_marker_genes: Whether to keep only marker genes in the expression matrix.
    :type keep_only_marker_genes: bool, False by default
    :param save_intermediate_object: Whether to save the intermediate versions of the Anndata object to the output directory.
    Allows for the script to be restarted from a later step if previous run fails. Useful when creating configs for large datasets.
    :type: bool
    """
    print(f"Converting CellBrowser config for project {project_name} to Anndata-Zarr object")
    config_converter = CellBrowserToAnndataZarrConverter(project_name, keep_only_marker_genes)
    config_is_valid = config_converter.download_config()
    if config_is_valid:
        config_converter.create_anndata_object()
        return config_converter.export_anndata_object()
    else:
        raise ValueError("CellBrowser config is not valid. Please check the error message above.")
