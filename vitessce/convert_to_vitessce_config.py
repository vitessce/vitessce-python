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

  def __init__(self, project, output_path):
     self.project = project
     self.output_path = output_path
     self.data = {}
     self.url_prefix = f"https://cells.ucsc.edu/{self.project}"
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
    if "topMarkers" in self.data and len(self.data["topMarkers"].keys()) > 0:
      print("Marker genes found. Filtering out all non-marker genes from Anndata object ...")
      marker_genes = [gene for sublist in self.data["topMarkers"].values() for gene in sublist]
      adata = adata[:, adata.var_names.isin(marker_genes)]

  def _write_to_zarr_store(self):
    data_dir = join('config_files', self.project)
    zarr_filepath = join(data_dir, "out.adata.zarr")

    obs_columns_list = adata.obs.columns.tolist()
    obsm_keys = list(adata.obsm.keys())
    var_cols_list = list(adata.var.keys())
    print(f"Obs columns: {obs_columns_list}")
    print(f"Obsm keys: {obsm_keys}")
    print(f"Var columns: {var_cols_list}")

    # This is done because of the optimize_adata function, which expects the obsm key to be numpy type
    for key in self.adata.obsm.keys():
        if isinstance(self.adata.obsm[key], pd.DataFrame):
            print(f"obsm {key} is an instance of DataFrame, converting it to numpy array.")
            self.adata.obsm[key] = adata.obsm[key].to_numpy()

    adata = optimize_adata(
        adata,
        obs_cols= obs_columns_list,
        obsm_keys=obsm_keys,
        optimize_X=True,
        var_cols=var_cols_list,
    )
    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])

  def _load_expr_matrix(self):
    gzip_file = None
    expr_matrix = None

    try:
      print("Downloading expression matrix ...")
      # Collect the exprMatrix and meta url suffixes
      expr_matrix_url_suffix = self.data["fileVersions"]["outMatrix"]["fname"].split("/")[-1]
    #   full_url = "/".join([self.url_prefix, expr_matrix_url_suffix])
    #   response = requests.get(full_url)
    #   response.raise_for_status()

    #   # Create a BytesIO object from the content
    #   gzip_file = io.BytesIO(response.content)
    #   print(f"Successfully downloaded expression matrix {full_url}.")
    # except Exception as e:
    #     print(f"Could not download expression matrix {full_url} because: {e}.")
    #     raise e
      full_path = os.path.join("config_files", self.project, expr_matrix_url_suffix)

      # Load the file
      with open(full_path, 'rb') as file:
          content = file.read()

      # If your file is gzipped, you can convert the content to a BytesIO object like this:
      gzip_file = io.BytesIO(content)
    except Exception as e:
      print(f"Could not download expression matrix because: {e}.")
      raise e

    with gzip.open(gzip_file, 'rt') as f:
      expr_matrix = pd.read_csv(f, sep='\t', index_col=0).T # transpose it, because of Scanpy
    
    # Now create anndata object
    self.adata = anndata.AnnData(X=expr_matrix)
    # TODO: hardcoding of gene
    self.adata.var['gene'] = self.adata.var_names
    first_gene = list(self.adata.var_names)[0]
    if len(first_gene.split("|")) == 2:
      # TODO: is there a way to have both?
      print("This dataset uses the format identifier|symbol for the ad.obs gene names (e.g. “ENSG0123123.3|HOX3”). We are keeping only the symbol.")
      self.adata.var_names = [x.split("|")[1] for x in list(self.adata.var_names)]

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

    for obj in self.data["coords"]:
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
        print("Success.")
      except Exception as e:
        print(f"Could not add {coord_type} to Anndata object because: {e}.")
        raise e
      print("Done adding coordinates to the Anndata object.")

  def _load_cell_metadata(self):
    try:
      print("Adding cell metadata to Anndata object ...")
      # # Load meta
      meta_url_suffix = self.data["fileVersions"]["outMeta"]["fname"].split("/")[-1]
      # response = requests.get("/".join([self.url_prefix, meta_url_suffix]))
      # response.raise_for_status()

      # # Create a BytesIO object from the content
      # meta_file = io.BytesIO(response.content)

      full_path = os.path.join("config_files", self.project, meta_url_suffix)

      # Load the file
      with open(full_path, 'rb') as file:
          content = file.read()

      # If your file is gzipped, you can convert the content to a BytesIO object like this:
      meta_file = io.BytesIO(content)

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
      "type" : "object",
      "properties" : {
          "fileVersions" : {
              "type" : "object",
              "properties" : {
                  "outMatrix" : {
                      "type" : "object",
                      "properties" : {
                          "fname" : {"type" : "string"},
                      },
                      "required": ["fname"]
                  },
              },
              "required": ["outMatrix"]
          },
          "coords" : {
              "type" : "array",
              "items": {
                  "type" : "object",
                  "properties" : {
                      "shortLabel" : {"type" : "string"},
                      "textFname" : {"type" : "string"},
                  },
                  "required": ["shortLabel", "textFname"]
              },
              "minItems": 1
          },
      },
      "required": ["fileVersions", "coords"]
    }

    try:
      validate(instance=self.data, schema=schema)
      print("JSON data is valid")
      return True
    except ValidationError as e:
      print("Invalid JSON data", e.message)
      return False

  def _download_config(self):  
    config_url = "/".join([self.url_prefix, "dataset.json"])
    try:
        # response = requests.get(config_url)
        # response.raise_for_status()
        # self.data = response.json()
        # todo: remove this when wanting to fetch the real data
        self.data = {
          "fileVersions": {
            "inMatrix": {
              "fname": "/hive/data/inside/cells/datasets/cortex-dev/exprMatrix.tsv.gz",
              "md5": "ad1f50ac39",
              "size": 124202312,
              "mtime": "2018-07-20 10:29:04"
            },
            "outMatrix": {
              "fname": "/usr/local/apache/htdocs-cells/cortex-dev/exprMatrix.tsv.gz",
              "md5": "645fd8e136",
              "size": 79656589,
              "mtime": "2021-09-13 14:50:29"
            },
            "inMeta": {
              "fname": "/hive/data/inside/cells/datasets/cortex-dev/meta.tsv",
              "md5": "a5649b68d8",
              "size": 181588,
              "mtime": "2018-09-06 17:00:46"
            },
            "outMeta": {
              "fname": "/usr/local/apache/htdocs-cells/cortex-dev/meta.tsv",
              "md5": "6e5a8aab3b",
              "size": 177326,
              "mtime": "2022-06-13 12:28:47"
            },
            "desc": {
              "fname": "desc.conf",
              "md5": "fdd764a0b3",
              "size": 1574,
              "mtime": "2020-06-08 15:33:49"
            },
            "conf": {
              "fname": "/hive/data/inside/cells/datasets/cortex-dev/cellbrowser.conf",
              "md5": "f7486b2d3c",
              "size": 2917,
              "mtime": "2022-06-13 12:28:42"
            }
          },
          "matrixArrType": "Float32",
          "sampleCount": 4261,
          "matrixWasFiltered": True,
          "coords": [
            {
              "name": "coords_0",
              "shortLabel": "t-SNE on WGCNA",
              "md5": "316cd1673f",
              "minX": 0,
              "maxX": 65535,
              "minY": 0,
              "maxY": 65535,
              "type": "Uint16",
              "textFname": "tMinusSNE_on_WGCNA.coords.tsv.gz",
              "labelMd5": "d41d8cd98f"
            },
            {
              "name": "coords_1",
              "shortLabel": "TriMap",
              "md5": "b3d8862cf3",
              "minX": 0,
              "maxX": 65535,
              "minY": 0,
              "maxY": 65535,
              "type": "Uint16",
              "textFname": "TriMap.coords.tsv.gz",
              "labelMd5": "d41d8cd98f"
            },
            {
              "name": "coords_2",
              "shortLabel": "T-SNE (scanpy)",
              "md5": "021a60ecb4",
              "minX": 0,
              "maxX": 65535,
              "minY": 1,
              "maxY": 65535,
              "type": "Uint16",
              "textFname": "TMinusSNE_scanpy.coords.tsv.gz",
              "labelMd5": "d41d8cd98f"
            },
            {
              "name": "coords_3",
              "shortLabel": "UMAP (scanpy)",
              "md5": "9d6dcf0ce6",
              "minX": 0,
              "maxX": 65535,
              "minY": 0,
              "maxY": 65535,
              "type": "Uint16",
              "textFname": "UMAP_scanpy.coords.tsv.gz",
              "labelMd5": "d41d8cd98f"
            },
            {
              "name": "coords_4",
              "shortLabel": "PAGA+ForceAtlas2 (scanpy)",
              "md5": "89d3a37436",
              "minX": 0,
              "maxX": 65534,
              "minY": 0,
              "maxY": 65535,
              "type": "Uint16",
              "textFname": "PAGAPlusForceAtlas2_scanpy.coords.tsv.gz",
              "labelMd5": "d41d8cd98f"
            }
          ],
        }
        print(f"Requested project name: {self.project}")
        print(f"Successfully fetched configuration: {config_url}.")
        return self._validate_config()
    except Exception as e:
        print(f"Could not get configuration for dataset {self.project} because: {e}.")
        return False

  def convert(self):
      print("Hello World!")
      config_is_valid = self._download_config()
      print(config_is_valid)
      if config_is_valid:
          self._load_expr_matrix()
          self._load_coordinates()
          self._load_cell_metadata()
          self._filter_data()
          self._write_to_zarr_store()