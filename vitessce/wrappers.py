import os
from os.path import join
import tempfile
import math

import numpy as np
import pandas as pd
import zarr
from numcodecs import Zlib
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from generate_tiff_offsets import get_offsets

from starlette.responses import JSONResponse, UJSONResponse
from starlette.routing import Route, Mount
from starlette.staticfiles import StaticFiles

from .constants import DataType as dt, FileType as ft
from .entities import Cells, CellSets, GenomicProfiles
from .routes import range_repsonse

class JsonRoute(Route):
    def __init__(self, path, endpoint, data_json):
        super().__init__(path, endpoint)
        self.data_json = data_json
    

class AbstractWrapper:
    """
    An abstract class that can be extended when
    implementing custom dataset object wrapper classes. 
    """

    def __init__(self, **kwargs):
        """
        Abstract constructor to be inherited by dataset wrapper classes.
        """
        pass

    def get_cells(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.CELLS` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_cell_sets(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.CELL_SETS` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_raster(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.RASTER` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_molecules(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.MOLECULES` data type.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_neighborhoods(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.NEIGHBORHOODS` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_expression_matrix(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.EXPRESSION_MATRIX` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_genomic_profiles(self, base_url, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the :class:`~vitessce.constants.DataType.GENOMIC_PROFILES` data type.
        Used internally by :class:`~vitessce.widget.VitessceWidget`.

        :param str base_url: The web server base url.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def _create_response_json(self, data_json):
        """
        Helper function that can be used for creating JSON responses.

        :param dict data_json: The data to return as JSON in the response body.
        :returns: The response handler function.
        :rtype: function
        """
        async def response_func(req):
            return UJSONResponse(data_json)
        return response_func

    def _get_data(self, data_type, base_url, dataset_uid, obj_i):
        if data_type == dt.CELLS:
            return self.get_cells(base_url, dataset_uid, obj_i)
        elif data_type == dt.CELL_SETS:
            return self.get_cell_sets(base_url, dataset_uid, obj_i)
        elif data_type == dt.RASTER:
            return self.get_raster(base_url, dataset_uid, obj_i)
        elif data_type == dt.MOLECULES:
            return self.get_molecules(base_url, dataset_uid, obj_i)
        elif data_type == dt.NEIGHBORHOODS:
            return self.get_neighborhoods(base_url, dataset_uid, obj_i)
        elif data_type == dt.EXPRESSION_MATRIX:
            return self.get_expression_matrix(base_url, dataset_uid, obj_i)
        elif data_type == dt.GENOMIC_PROFILES:
            return self.get_genomic_profiles(base_url, dataset_uid, obj_i)

    def _get_url(self, base_url, dataset_uid, obj_i, suffix):
        return base_url + self._get_route(dataset_uid, obj_i, suffix)

    def _get_route(self, dataset_uid, obj_i, suffix):
        return f"/{dataset_uid}/{obj_i}/{suffix}"

class MultiImageWrapper(AbstractWrapper):
    """
    Wrap multiple imaging datasets by creating an instance of the ``MultiImageWrapper`` class.

    :param list image_wrappers: A list of imaging wrapper classes (only :class:`~vitessce.wrappers.OmeTiffWrapper` supported now)
    :param \*\*kwargs: Keyword arguments inherited from :class:`~vitessce.wrappers.AbstractWrapper`
    """
    def __init__(self, image_wrappers, **kwargs):
        super().__init__(**kwargs)
        self.image_wrappers = image_wrappers

    def create_raster_json(self, base_url="", dataset_uid="", obj_i=""):
        raster_json = {
            "schemaVersion": "0.0.2",
            "images": [],
            "renderLayers": []
        }
        for image in self.image_wrappers:
            image_json = image.create_image_json(
                image.get_img_url(base_url, dataset_uid, obj_i),
                image.get_offsets_url(base_url, dataset_uid, obj_i)
            )
            raster_json['images'].append(image_json)
            raster_json['renderLayers'].append(image.name)
        return raster_json
    
    def get_raster(self, base_url="", dataset_uid="", obj_i=""):
        raster_json = self.create_raster_json(base_url, dataset_uid, obj_i)
        obj_routes = []
        for image in self.image_wrappers:
            obj_routes = obj_routes + image.get_routes(base_url, dataset_uid, obj_i)
        obj_file_defs = [
            {
                "type": dt.RASTER.value,
                "fileType": ft.RASTER_JSON.value,
                "options": raster_json
            }
        ]

        return obj_file_defs, obj_routes

class OmeTiffWrapper(AbstractWrapper):

    """
    Wrap an OME-TIFF File by creating an instance of the ``OmeTiffWrapper`` class.

    :param str img_path: A local filepath to an OME-TIFF file.
    :param str img_url: A remote URL of an OME-TIFF file.
    :param str name: The display name for this OME-TIFF within Vitessce.
    :param list[number] transformation_matrix: A column-major ordered matrix for transforming this image (see http://www.opengl-tutorial.org/beginners-tutorials/tutorial-3-matrices/#homogeneous-coordinates for more information).
    :param \*\*kwargs: Keyword arguments inherited from :class:`~vitessce.wrappers.AbstractWrapper`
    """

    def __init__(self, img_path="", img_url="", offsets_url="", name="", transformation_matrix=None, **kwargs):
        super().__init__(**kwargs)
        self.name = name
        self._img_path = img_path
        self._img_url = img_url
        self._offsets_url = offsets_url
        self._transformation_matrix = transformation_matrix

    def create_raster_json(self, img_url, offsets_url=""):
        raster_json = {
            "schemaVersion": "0.0.2",
            "images": [self.create_image_json(img_url, offsets_url)],
        }
        return raster_json
    
    def create_image_json(self, img_url, offsets_url=""):
        metadata = {}
        image = {
            "name": self.name,
            "type": "ome-tiff",
            "url": img_url,
        }
        if offsets_url != "":
            metadata["omeTiffOffsetsUrl"] = offsets_url
        if self._transformation_matrix is not None:
            metadata["transform"] = {
                "matrix": self._transformation_matrix
            }
        # Only attach metadata if there is some - otherwise schema validation fails.
        if len(metadata.keys()) > 0:
            image["metadata"] = metadata
        return image
    
    def _get_image_dir(self):
        return os.path.dirname(self._img_path)
    
    def _get_img_filename(self):
        return os.path.basename(self._img_path)

    def get_img_url(self, base_url="", dataset_uid="", obj_i=""):
        if self._img_url != "":
            return self._img_url
        img_url = self._get_url(base_url, dataset_uid, obj_i, self._get_img_filename())
        return img_url

    def get_offsets_path_name(self):
        return f"{self._get_img_filename().split('.')[0]}.offsets.json"
    
    def get_offsets_url(self, base_url="", dataset_uid="", obj_i=""):
        if self._offsets_url != "" or self._img_url != "":
            return self._offsets_url
        offsets_url = self._get_url(base_url, dataset_uid, obj_i, self.get_offsets_path_name())
        return offsets_url
    
    def get_routes(self, base_url="", dataset_uid="", obj_i=""):
        obj_routes = [
            Route(self._get_route(dataset_uid, obj_i, self._get_img_filename()), lambda req: range_repsonse(req, self._img_path))
        ]
        if self._img_path != "":
            offsets = get_offsets(self._img_path)
            obj_routes.append(
                JsonRoute(self._get_route(dataset_uid, obj_i, self.get_offsets_path_name()),
                        self._create_response_json(offsets), offsets),
            )
        return obj_routes

    def get_raster(self, base_url="", dataset_uid="", obj_i=""):
        obj_routes = self.get_routes(base_url, dataset_uid, obj_i)
        img_url = self.get_img_url(base_url, dataset_uid, obj_i)
        offsets_url = self.get_offsets_url(base_url, dataset_uid, obj_i)
        raster_json = self.create_raster_json(img_url, offsets_url)
        obj_file_defs = [
            {
                "type": dt.RASTER.value,
                "fileType": ft.RASTER_JSON.value,
                "options": raster_json
            }
        ]

        return obj_file_defs, obj_routes


class OmeZarrWrapper(AbstractWrapper):

    def __init__(self, z, name="", **kwargs):
        super().__init__(**kwargs)
        self.z = z
        self.name = name

    def create_raster_json(self, img_url):
        raster_json = {
            "schemaVersion": "0.0.2",
            "images": [
                {
                    "name": self.name,
                    "type": "zarr",
                    "url": img_url,
                    "metadata": {
                        "dimensions": [
                            {
                                "field": "channel",
                                "type": "nominal",
                                "values": [
                                    "DAPI - Hoechst (nuclei)",
                                    "FITC - Laminin (basement membrane)",
                                    "Cy3 - Synaptopodin (glomerular)",
                                    "Cy5 - THP (thick limb)"
                                ]
                            },
                            {
                                "field": "y",
                                "type": "quantitative",
                                "values": None
                            },
                            {
                                "field": "x",
                                "type": "quantitative",
                                "values": None
                            }
                        ],
                        "isPyramid": True,
                        "transform": {
                            "scale": 1,
                            "translate": {
                                "x": 0,
                                "y": 0,
                            }
                        }
                    }
                }
            ],
        }
        return raster_json

    def get_raster(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        if type(self.z) == zarr.hierarchy.Group:
            img_dir_path = self.z.store.path

            raster_json = self.create_raster_json(
                self._get_url(base_url, dataset_uid, obj_i, "raster_img"),
            )

            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "raster_img"),
                        app=StaticFiles(directory=img_dir_path, html=False)),
                JsonRoute(self._get_route(dataset_uid, obj_i, "raster"),
                        self._create_response_json(raster_json), raster_json)
            ]
            obj_file_defs = [
                {
                    "type": dt.RASTER.value,
                    "fileType": ft.RASTER_JSON.value,
                    "url": self._get_url(base_url, dataset_uid, obj_i, "raster")
                }
            ]

        return obj_file_defs, obj_routes


class AnnDataWrapper(AbstractWrapper):
    def __init__(self, adata, use_highly_variable_genes=True, cell_set_obs_cols=None, spatial_obsm_key=None, **kwargs):
        """
        Wrap an AnnData object by creating an instance of the ``AnnDataWrapper`` class.

        :param adata: An AnnData object containing single-cell experiment data.
        :type adata: anndata.AnnData
        :param bool use_highly_variable_genes: When creating outputs with genes, should only the genes marked as highly variable be used?
        :param list[str] cell_set_obs_cols: A list of column names of the ``adata.obs`` dataframe that should be used for creating cell sets.
        :param \*\*kwargs: Keyword arguments inherited from :class:`~vitessce.wrappers.AbstractWrapper`
        """
        super().__init__(**kwargs)
        self.adata = adata
        self.tempdir = tempfile.mkdtemp()
        self.use_highly_variable_genes = use_highly_variable_genes
        self.cell_set_obs_cols = cell_set_obs_cols
        self.spatial_obsm_key = spatial_obsm_key

    def create_cells_json(self):
        adata = self.adata
        available_embeddings = set(adata.obsm.keys()) - set([self.spatial_obsm_key])

        cell_ids = adata.obs.index.tolist()
        cells = Cells(cell_ids=cell_ids)
        for e in available_embeddings:
            mapping = adata.obsm[e][:, 0:2].tolist()
            cells.add_mapping(e, mapping)
        
        if self.spatial_obsm_key is not None:
            centroids = adata.obsm[self.spatial_obsm_key][:, 0:2].tolist()
            cells.add_centroids(centroids)
        
        return cells.json

    def create_cell_sets_json(self):
        adata = self.adata

        cell_set_obs_cols = self.cell_set_obs_cols

        cell_sets = CellSets()

        if cell_set_obs_cols is not None and len(cell_set_obs_cols) > 0:
            # Each `cell_set_obs_col` is a column name in the `adata.obs` dataframe,
            # which we want to turn into a hierarchy of cell sets.
            for cell_set_obs_col in cell_set_obs_cols:
                cell_sets.add_level_zero_node(cell_set_obs_col)

                cell_ids = adata.obs.index.tolist()
                cluster_ids = adata.obs[cell_set_obs_col].unique().tolist()
                cell_cluster_ids = adata.obs[cell_set_obs_col].values.tolist()

                cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))

                for cluster_id in sorted(cluster_ids):
                    cell_set = [
                        str(cell_id)
                        for cell_id, cell_cluster_id in cell_cluster_tuples
                        if cell_cluster_id == cluster_id
                    ]
                    cell_sets.add_node(str(cluster_id), [cell_set_obs_col], cell_set)
            return cell_sets.json
        return None
    
    def create_exp_matrix_zarr(self, zarr_filepath):
        adata = self.adata
        gexp_arr = adata.X

        cell_list = adata.obs.index.values.tolist()
        gene_list = adata.var.index.values.tolist()

        if type(gexp_arr) == csr_matrix:
            # Convert from SciPy sparse format to NumPy dense format
            gexp_arr = gexp_arr.toarray()
        
        if self.use_highly_variable_genes and 'highly_variable' in adata.var.columns.values.tolist():
            # Restrict the gene expression matrix to only the genes marked as highly variable
            gene_list = adata.var.index[adata.var['highly_variable']].values.tolist()
            gexp_arr = gexp_arr[:,adata.var['highly_variable'].values]

        
        # Re-scale the gene expression values between 0 and 255
        gexp_arr_min = gexp_arr.min()
        gexp_arr_max = gexp_arr.max()
        gexp_arr_range = gexp_arr_max - gexp_arr_min
        gexp_arr_ratio = 255 / gexp_arr_range

        gexp_norm_arr = (gexp_arr - gexp_arr_min) * gexp_arr_ratio
    
        z = zarr.open(
            zarr_filepath,
            mode='w',
            shape=gexp_norm_arr.shape,
            dtype='uint8',
            compressor=Zlib(level=1)
        )

        z[:] = gexp_norm_arr
        # observations: cells (rows)
        z.attrs["rows"] = cell_list
        # variables: genes (columns)
        z.attrs["cols"] = gene_list
        
        return

    def get_cells(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        cells_json = self.create_cells_json()

        obj_routes = [
            JsonRoute(self._get_route(dataset_uid, obj_i, "cells"),
                    self._create_response_json(cells_json), cells_json),
        ]
        obj_file_defs = [
            {
                "type": dt.CELLS.value,
                "fileType": ft.CELLS_JSON.value,
                "url": self._get_url(base_url, dataset_uid, obj_i, "cells")
            }
        ]

        return obj_file_defs, obj_routes

    def get_cell_sets(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
  
        cell_sets_json = self.create_cell_sets_json()

        if cell_sets_json is not None:
            obj_routes = [
                JsonRoute(self._get_route(dataset_uid, obj_i, "cell-sets"),
                        self._create_response_json(cell_sets_json), cell_sets_json),
            ]
            obj_file_defs = [
                {
                    "type": dt.CELL_SETS.value,
                    "fileType": ft.CELL_SETS_JSON.value,
                    "url": self._get_url(base_url, dataset_uid, obj_i, "cell-sets")
                }
            ]

        return obj_file_defs, obj_routes
    
    def get_expression_matrix(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        zarr_tempdir = self.tempdir
        zarr_filepath = join(zarr_tempdir, 'matrix.zarr')

        self.create_exp_matrix_zarr(zarr_filepath)

        if zarr_tempdir is not None:
            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "expression"),
                    app=StaticFiles(directory=os.path.dirname(zarr_filepath), html=False, check_dir=False)),
            ]

            obj_file_defs = [
                {
                    "type": dt.EXPRESSION_MATRIX.value,
                    "fileType": ft.EXPRESSION_MATRIX_ZARR.value,
                    "url": self._get_url(base_url, dataset_uid, obj_i, "expression/matrix.zarr"),
                }
            ]

        return obj_file_defs, obj_routes

class SnapWrapper(AbstractWrapper):

    # The Snap file is difficult to work with.
    # For now we can use the processed cell-by-bin MTX file
    # However, the HuBMAP pipeline currently computes this with resolution 5000
    # TODO: Make a PR to sc-atac-seq-pipeline to output this at a higher resolution (e.g. 200)
    # https://github.com/hubmapconsortium/sc-atac-seq-pipeline/blob/develop/bin/snapAnalysis.R#L93

    def __init__(self, in_mtx, in_barcodes_df, in_bins_df, in_clusters_df, starting_resolution=5000, **kwargs):
        super().__init__(**kwargs)
        self.in_mtx = in_mtx # scipy.sparse.coo.coo_matrix (filtered_cell_by_bin.mtx)
        self.in_barcodes_df = in_barcodes_df # pandas dataframe (barcodes.txt)
        self.in_bins_df = in_bins_df # pandas dataframe (bins.txt)
        self.in_clusters_df = in_clusters_df # pandas dataframe (umap_coords_clusters.csv)

        self.tempdir = tempfile.mkdtemp()

        self.starting_resolution = starting_resolution

        # Convert to dense matrix if sparse.
        if type(in_mtx) == coo_matrix:
            self.in_mtx = in_mtx.toarray()


    def create_genomic_multivec_zarr(self, zarr_filepath):
        in_mtx = self.in_mtx
        in_clusters_df = self.in_clusters_df
        in_barcodes_df = self.in_barcodes_df
        in_bins_df = self.in_bins_df

        starting_resolution = self.starting_resolution

        # The bin datafram consists of one column like chrName:binStart-binEnd
        def convert_bin_name_to_chr_name(bin_name):
            try:
                return bin_name[:bin_name.index(':')]
            except ValueError:
                return np.nan
        def convert_bin_name_to_chr_start(bin_name):
            try:
                return int(bin_name[bin_name.index(':')+1:bin_name.index('-')])
            except ValueError:
                return np.nan
        def convert_bin_name_to_chr_end(bin_name):
            try:
                return int(bin_name[bin_name.index('-')+1:])
            except ValueError:
                return np.nan
        
        # The genome assembly is GRCh38 but the chromosome names in the bin names do not start with the "chr" prefix.
        # This is incompatible with the chromosome names from `negspy`, so we need to append the prefix.
        in_bins_df[0] = in_bins_df[0].apply(lambda x: "chr" + x)
        
        in_bins_df["chr_name"] = in_bins_df[0].apply(convert_bin_name_to_chr_name)
        in_bins_df["chr_start"] = in_bins_df[0].apply(convert_bin_name_to_chr_start)
        in_bins_df["chr_end"] = in_bins_df[0].apply(convert_bin_name_to_chr_end)

        # Drop any rows that had incorrect bin strings (missing a chromosome name, bin start, or bin end value).
        in_bins_df = in_bins_df.dropna(subset=["chr_name", "chr_start", "chr_end"]).copy()

        # Ensure that the columns have the expect types.
        in_bins_df["chr_name"] = in_bins_df["chr_name"].astype(str)
        in_bins_df["chr_start"] = in_bins_df["chr_start"].astype(int)
        in_bins_df["chr_end"] = in_bins_df["chr_end"].astype(int)

        # Create the Zarr store for the outputs.
        out_f = zarr.open(zarr_filepath, mode='w')

        # Get a list of clusters.
        in_clusters_df["cluster"] = in_clusters_df["cluster"].astype(str)
        cluster_ids = in_clusters_df["cluster"].unique().tolist()
        cluster_ids.sort(key=int)

        cluster_paths = [ [ "Clusters", cluster_id ] for cluster_id in cluster_ids ]
        
        # "SnapTools performs quantification using a specified aligner, and HuBMAP has standardized on BWA with the GRCh38 reference genome"
        # Reference: https://github.com/hubmapconsortium/sc-atac-seq-pipeline/blob/bb023f95ca3330128bfef41cc719ffcb2ee6a190/README.md
        genomic_profiles = GenomicProfiles(out_f, profile_paths=cluster_paths, assembly='hg38', starting_resolution=starting_resolution)
        chrom_name_to_length = genomic_profiles.chrom_name_to_length

        # Create each chromosome dataset.
        for chr_name, chr_len in chrom_name_to_length.items():
            # The bins dataframe frustratingly does not contain every bin.
            # We need to figure out which bins are missing.

            # We want to check for missing bins in each chromosome separately,
            # otherwise too much memory is used during the join step.
            chr_bins_in_df = in_bins_df.loc[in_bins_df["chr_name"] == chr_name]
            if chr_bins_in_df.shape[0] == 0:
                # No processing or output is necessary if there is no data for this chromosome.
                # Continue on through all resolutions of this chromosome to the next chromosome.
                continue
            
            # Determine the indices of the matrix at which the bins for this chromosome start and end.
            chr_bin_i_start = int(chr_bins_in_df.head(1).iloc[0].name)
            chr_bin_i_end = int(chr_bins_in_df.tail(1).iloc[0].name) + 1
            
            # Extract the part of the matrix corresponding to the current chromosome.
            chr_mtx = in_mtx[:,chr_bin_i_start:chr_bin_i_end]

            # Create a list of the "ground truth" bins (all bins from position 0 to the end of the chromosome).
            # We will join the input bins onto this dataframe to determine which bins are missing.
            chr_bins_gt_df = pd.DataFrame()
            chr_bins_gt_df["chr_start"] = np.arange(0, math.ceil(chr_len/starting_resolution)) * starting_resolution
            chr_bins_gt_df["chr_end"] = chr_bins_gt_df["chr_start"] + starting_resolution
            chr_bins_gt_df["chr_start"] = chr_bins_gt_df["chr_start"] + 1
            chr_bins_gt_df["chr_start"] = chr_bins_gt_df["chr_start"].astype(int)
            chr_bins_gt_df["chr_end"] = chr_bins_gt_df["chr_end"].astype(int)
            chr_bins_gt_df["chr_name"] = chr_name
            chr_bins_gt_df[0] = chr_bins_gt_df.apply(lambda r: f"{r['chr_name']}:{r['chr_start']}-{r['chr_end']}", axis='columns')
            
            # We will add a new column "i", which should match the _old_ index, so that we will be able join with the data matrix on the original indices.
            # For the new rows, we will add values for the "i" column that are greater than any of the original indices,
            # to prevent any joining with the incoming data matrix onto these bins for which the data is missing.
            chr_bins_in_df = chr_bins_in_df.reset_index(drop=True)
            chr_bins_in_df["i"] = chr_bins_in_df.index.values
            chr_bins_gt_df["i"] = chr_bins_gt_df.index.values + (in_mtx.shape[1] + 1)
            
            # Set the full bin string column as the index of both data frames.
            chr_bins_gt_df = chr_bins_gt_df.set_index(0)
            chr_bins_in_df = chr_bins_in_df.set_index(0)
            
            # Join the input bin subset dataframe right onto the full bin ground truth dataframe.
            chr_bins_in_join_df = chr_bins_in_df.join(chr_bins_gt_df, how='right', lsuffix="", rsuffix="_gt")
            # The bins which were not present in the input will have NaN values in the "i" column.
            # For these rows, we replace the NaN values with the much higher "i_gt" values which will not match to any index of the data matrix.
            chr_bins_in_join_df["i"] = chr_bins_in_join_df.apply(lambda r: r['i'] if pd.notna(r['i']) else r['i_gt'], axis='columns').astype(int)

            # Clean up the joined data frame by removing unnecessary columns.
            chr_bins_in_join_df = chr_bins_in_join_df.drop(columns=['chr_name', 'chr_start', 'chr_end', 'i_gt'])
            chr_bins_in_join_df = chr_bins_in_join_df.rename(columns={'chr_name_gt': 'chr_name', 'chr_start_gt': 'chr_start', 'chr_end_gt': 'chr_end'})
            
            # Create a dataframe from the data matrix, so that we can join to the joined bins dataframe.
            chr_mtx_df = pd.DataFrame(data=chr_mtx.T)
            
            chr_bins_i_df = chr_bins_in_join_df.drop(columns=['chr_name', 'chr_start', 'chr_end'])

            # Join the data matrix dataframe and the bins dataframe.
            # Bins that are missing from the data matrix will have "i" values higher than any of the data matrix dataframe row indices,
            # and therefore the data values for these bins in the resulting joined dataframe will all be NaN.
            chr_mtx_join_df = chr_bins_i_df.join(chr_mtx_df, how='left', on='i')
            # We fill in these NaN values with 0.
            chr_mtx_join_df = chr_mtx_join_df.fillna(value=0.0)
            
            # Drop the "i" column, since it is not necessary now that we have done the join.
            chr_mtx_join_df = chr_mtx_join_df.drop(columns=['i'])
            # Obtain the new full data matrix, which contains values for all bins of the chromosome.
            chr_mtx = chr_mtx_join_df.values.T

            # Fill in the Zarr store with data for each cluster.
            for cluster_index, cluster_id in enumerate(cluster_ids):
                # Get the list of cells in the current cluster.
                cluster_df = in_clusters_df.loc[in_clusters_df["cluster"] == cluster_id]
                cluster_cell_ids = cluster_df.index.values.tolist()
                cluster_num_cells = len(cluster_cell_ids)
                cluster_cells_tf = (in_barcodes_df[0].isin(cluster_cell_ids)).values

                # Get the rows of the data matrix corresponding to the cells in this cluster.
                cluster_cell_by_bin_mtx = chr_mtx[cluster_cells_tf,:]
                # Take the sum of this cluster along the cells axis.
                cluster_profile = cluster_cell_by_bin_mtx.sum(axis=0)

                genomic_profiles.add_profile(cluster_profile, chr_name, cluster_index)
        
        return

    def get_genomic_profiles(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
        
        zarr_tempdir = self.tempdir
        zarr_filepath = join(zarr_tempdir, 'profiles.zarr')

        print("Please wait, the following conversion is slow")
        self.create_genomic_multivec_zarr(zarr_filepath)

        if zarr_tempdir is not None:
            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "genomic"),
                    app=StaticFiles(directory=os.path.dirname(zarr_filepath), html=False, check_dir=False)),
            ]

            obj_file_defs = [
                {
                    "type": dt.GENOMIC_PROFILES.value,
                    "fileType": ft.GENOMIC_PROFILES_ZARR.value,
                    "url": self._get_url(base_url, dataset_uid, obj_i, "genomic/profiles.zarr")
                }
            ]

        return obj_file_defs, obj_routes
    

    def create_cell_sets_json(self):
        in_clusters_df = self.in_clusters_df
        cell_sets = CellSets()
        cell_sets.add_level_zero_node('Clusters')

        cell_ids = in_clusters_df.index.values.tolist()
        in_clusters_df['cluster'] = in_clusters_df['cluster'].astype(str)
        cluster_ids = in_clusters_df['cluster'].unique().tolist()
        cluster_ids.sort(key=int)
        cell_cluster_ids = in_clusters_df['cluster'].values.tolist()

        cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))

        for cluster_id in cluster_ids:
            cell_set = [
                str(cell_id)
                for cell_id, cell_cluster_id in cell_cluster_tuples
                if cell_cluster_id == cluster_id
            ]
            cell_sets.add_node(str(cluster_id), ['Clusters'], cell_set)

        return cell_sets.json
    
    def get_cell_sets(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        cell_sets_json = self.create_cell_sets_json()

        obj_routes = [
            JsonRoute(self._get_route(dataset_uid, obj_i, "cell-sets"),
                    self._create_response_json(cell_sets_json), cell_sets_json),
        ]
        obj_file_defs = [
            {
                "type": dt.CELL_SETS.value,
                "fileType": ft.CELL_SETS_JSON.value,
                "url": self._get_url(base_url, dataset_uid, obj_i, "cell-sets")
            }
        ]

        return obj_file_defs, obj_routes
    
    def create_cells_json(self):
        in_clusters_df = self.in_clusters_df

        cell_ids = in_clusters_df.index.tolist()
        cells = Cells(cell_ids=cell_ids)
        mapping = in_clusters_df[["umap.1", "umap.2"]].values.tolist()
        cells.add_mapping("UMAP", mapping)
        return cells.json
    
    def get_cells(self, base_url, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        cells_json = self.create_cells_json()

        obj_routes = [
            JsonRoute(self._get_route(dataset_uid, obj_i, "cells"),
                    self._create_response_json(cells_json), cells_json),
        ]
        obj_file_defs = [
            {
                "type": dt.CELLS.value,
                "fileType": ft.CELLS_JSON.value,
                "url": self._get_url(base_url, dataset_uid, obj_i, "cells")
            }
        ]

        return obj_file_defs, obj_routes
