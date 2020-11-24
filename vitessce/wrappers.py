import os
from os.path import join
import tempfile
import math

import numpy as np
import zarr
from numcodecs import Zlib
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

from starlette.responses import JSONResponse, UJSONResponse
from starlette.routing import Route, Mount
from starlette.staticfiles import StaticFiles

from .constants import DataType as dt, FileType as ft

class AbstractWrapper:
    """
    An abstract class that can be extended when
    implementing custom dataset object wrapper classes. 
    """
    def get_cells(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``cells`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_cell_sets(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``cell-sets`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_raster(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``raster`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_molecules(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``molecules`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_neighborhoods(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``neighborhoods`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_expression_matrix(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``expression-matrix`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
        :param str dataset_uid: The unique identifier for the dataset parent of this data object.
        :param int obj_i: The index of this data object child within its dataset parent.

        :returns: The file definitions and server routes.
        :rtype: tuple[list[dict], list[starlette.routing.Route]]
        """
        raise NotImplementedError()

    def get_genomic_profiles(self, port, dataset_uid, obj_i):
        """
        Get the file definitions and server routes
        corresponding to the ``genomic-profiles`` data type.

        :param int port: The web server port, meant to be used in the localhost URLs in the file definitions.
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

    def _get_data(self, data_type, port, dataset_uid, obj_i):
        if data_type == dt.CELLS:
            return self.get_cells(port, dataset_uid, obj_i)
        elif data_type == dt.CELL_SETS:
            return self.get_cell_sets(port, dataset_uid, obj_i)
        elif data_type == dt.RASTER:
            return self.get_raster(port, dataset_uid, obj_i)
        elif data_type == dt.MOLECULES:
            return self.get_molecules(port, dataset_uid, obj_i)
        elif data_type == dt.NEIGHBORHOODS:
            return self.get_neighborhoods(port, dataset_uid, obj_i)
        elif data_type == dt.EXPRESSION_MATRIX:
            return self.get_expression_matrix(port, dataset_uid, obj_i)
        elif data_type == dt.GENOMIC_PROFILES:
            return self.get_genomic_profiles(port, dataset_uid, obj_i)

    def _get_url(self, port, dataset_uid, obj_i, suffix):
        return f"http://localhost:{port}/{dataset_uid}/{obj_i}/{suffix}"

    def _get_route(self, dataset_uid, obj_i, suffix):
        return f"/{dataset_uid}/{obj_i}/{suffix}"


class OmeTiffWrapper(AbstractWrapper):

    def __init__(self, img_path, offsets_path=None, name=""):
        self.img_path = img_path
        self.offsets_path = offsets_path
        self.name = name

    def _create_raster_json(self, img_url, offsets_url):
        raster_json = {
            "schemaVersion": "0.0.2",
            "images": [
                {
                    "name": self.name,
                    "type": "ome-tiff",
                    "url": img_url,
                    "metadata": {
                        **({
                            "omeTiffOffsetsUrl": offsets_url,
                        } if offsets_url is not None else {})
                    }
                }
            ],
        }
        return raster_json

    def _get_offsets_dir(self):
        return os.path.dirname(self.offsets_path)
    
    def _get_offsets_filename(self):
        return os.path.basename(self.offsets_path)

    def get_raster(self, port, dataset_uid, obj_i):
        img_dir_path, img_url = self.img_path, self._get_url(port, dataset_uid, obj_i, "raster_img")
        offsets_dir_path, offsets_url = (None, None) if self.offsets_path is None else (self._get_offsets_dir(), self._get_url(port, dataset_uid, obj_i, join("raster_offsets", self._get_offsets_filename())))

        raster_json = self._create_raster_json(img_url, offsets_url)

        obj_routes = [
            Mount(self._get_route(dataset_uid, obj_i, "raster_img"),
                  app=StaticFiles(directory=img_dir_path, html=False, check_dir=False)),
            Route(self._get_route(dataset_uid, obj_i, "raster"),
                  self._create_response_json(raster_json))
        ]
        if self.offsets_path is not None:
            obj_routes.append(
                Mount(self._get_route(dataset_uid, obj_i, "raster_offsets"),
                      app=StaticFiles(directory=offsets_dir_path, html=False, check_dir=False))
            )

        obj_file_defs = [
            {
                "type": dt.RASTER.value,
                "fileType": ft.RASTER_JSON.value,
                "url": self._get_url(port, dataset_uid, obj_i, "raster")
            }
        ]

        return obj_file_defs, obj_routes


class ZarrDirectoryStoreWrapper(AbstractWrapper):

    def __init__(self, z, name=""):
        self.z = z
        self.name = name

    def _create_raster_json(self, img_url):
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

    def get_raster(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        if type(self.z) == zarr.hierarchy.Group:
            img_dir_path = self.z.store.path

            raster_json = self._create_raster_json(
                self._get_url(port, dataset_uid, obj_i, "raster_img"),
            )

            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "raster_img"),
                        app=StaticFiles(directory=img_dir_path, html=False)),
                Route(self._get_route(dataset_uid, obj_i, "raster"),
                        self._create_response_json(raster_json))
            ]
            obj_file_defs = [
                {
                    "type": dt.RASTER.value,
                    "fileType": ft.RASTER_JSON.value,
                    "url": self._get_url(port, dataset_uid, obj_i, "raster")
                }
            ]

        return obj_file_defs, obj_routes


class AnnDataWrapper(AbstractWrapper):
    def __init__(self, adata, use_highly_variable_genes=True):
        self.adata = adata
        self.tempdir = tempfile.mkdtemp()

        self.use_highly_variable_genes = use_highly_variable_genes

    def _create_cells_json(self):
        adata = self.adata
        available_embeddings = list(adata.obsm.keys())

        cell_ids = adata.obs.index.tolist()
        cell_mappings = []
        for e in available_embeddings:
            mapping = adata.obsm[e][:, 0:2].tolist()
            cell_mappings.append(list(zip(
                [e for i in range(len(mapping))],
                mapping
            )))
        cell_mappings_zip = list(zip(*cell_mappings))
        cells_json = dict(zip(
            cell_ids,
            [
                {'mappings': dict(cell_mapping), 'genes': {}}
                for cell_mapping in cell_mappings_zip
            ]
        ))
        return cells_json

    def _create_cell_sets_json(self):
        adata = self.adata
        cell_sets_json = {
            "datatype": "cell",
            "version": "0.1.2",
            "tree": [{
                "name": "Clusters",
                "children": []
            }]
        }

        cell_ids = adata.obs.index.tolist()
        cluster_ids = adata.obs['CellType'].unique().tolist()
        cell_cluster_ids = adata.obs['CellType'].values.tolist()

        cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))

        for cluster_id in cluster_ids:
            cell_sets_json["tree"][0]["children"].append({
                "name": str(cluster_id),
                "set": [
                    str(cell_id)
                    for cell_id, cell_cluster_id in cell_cluster_tuples
                    if cell_cluster_id == cluster_id
                ]
            })

        return cell_sets_json
    
    def _create_exp_matrix_zarr(self, zarr_filepath):
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

    def get_cells(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        cells_json = self._create_cells_json()

        obj_routes = [
            Route(self._get_route(dataset_uid, obj_i, "cells"),
                    self._create_response_json(cells_json)),
        ]
        obj_file_defs = [
            {
                "type": dt.CELLS.value,
                "fileType": ft.CELLS_JSON.value,
                "url": self._get_url(port, dataset_uid, obj_i, "cells")
            }
        ]

        return obj_file_defs, obj_routes

    def get_cell_sets(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

            
        cell_sets_json = self._create_cell_sets_json()

        obj_routes = [
            Route(self._get_route(dataset_uid, obj_i, "cell-sets"),
                    self._create_response_json(cell_sets_json)),
        ]
        obj_file_defs = [
            {
                "type": dt.CELL_SETS.value,
                "fileType": ft.CELL_SETS_JSON.value,
                "url": self._get_url(port, dataset_uid, obj_i, "cell-sets")
            }
        ]

        return obj_file_defs, obj_routes
    
    def get_expression_matrix(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        zarr_tempdir = self.tempdir
        zarr_filepath = join(zarr_tempdir, 'matrix.zarr')

        self._create_exp_matrix_zarr(zarr_filepath)

        if zarr_tempdir is not None:
            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "expression"),
                    app=StaticFiles(directory=os.path.dirname(zarr_filepath), html=False, check_dir=False)),
            ]

            obj_file_defs = [
                {
                    "type": dt.EXPRESSION_MATRIX.value,
                    "fileType": ft.EXPRESSION_MATRIX_ZARR.value,
                    "url": self._get_url(port, dataset_uid, obj_i, "expression/matrix.zarr")
                }
            ]

        return obj_file_defs, obj_routes
        


class LoomWrapper(AbstractWrapper):

    def __init__(self, loom):
        self.loom = loom

    def get_cells(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
        """
        try:
            import loompy
            if type(self.loom) == loompy.LoomConnection:
                pass
                # TODO: append routes
                # TODO: add file definitions
        except ImportError:
            pass
        """
        return obj_file_defs, obj_routes

class SnapToolsWrapper(AbstractWrapper):

    # The Snap file is difficult to work with.
    # For now we can use the processed cell-by-bin MTX file
    # However, the HuBMAP pipeline currently computes this with resolution 5000
    # https://github.com/hubmapconsortium/sc-atac-seq-pipeline/blob/develop/bin/snapAnalysis.R#L93

    def __init__(self, in_mtx, in_barcodes_df, in_bins_df, in_clusters_df, starting_resolution=5000):
        self.in_mtx = in_mtx # scipy.sparse.coo.coo_matrix (filtered_cell_by_bin.mtx)
        self.in_barcodes_df = in_barcodes_df # pandas dataframe (barcodes.txt)
        self.in_bins_df = in_bins_df # pandas dataframe (bins.txt)
        self.in_clusters_df = in_clusters_df # pandas dataframe (umap_coords_clusters.csv)

        self.tempdir = tempfile.mkdtemp()

        self.starting_resolution = starting_resolution

        # Convert to dense if sparse
        if type(in_mtx) == coo_matrix:
            self.in_mtx = in_mtx.toarray()


    def _create_genomic_multivec_zarr(self, zarr_filepath):
        in_mtx = self.in_mtx
        in_clusters_df = self.in_clusters_df
        in_barcodes_df = self.in_barcodes_df
        in_bins_df = self.in_bins_df

        starting_resolution = self.starting_resolution

        def convert_bin_name_to_chr_name(bin_name):
            return bin_name[:bin_name.index(':')]
        def convert_bin_name_to_chr_start(bin_name):
            return int(bin_name[bin_name.index(':')+1:bin_name.index('-')])
        def convert_bin_name_to_chr_end(bin_name):
            return int(bin_name[bin_name.index('-')+1:])
        
        in_bins_df["chr_name"] = in_bins_df[0].apply(convert_bin_name_to_chr_name).astype(str)
        in_bins_df["chr_start"] = in_bins_df[0].apply(convert_bin_name_to_chr_start).astype(int)
        in_bins_df["chr_end"] = in_bins_df[0].apply(convert_bin_name_to_chr_end).astype(int)

        out_f = zarr.open(zarr_filepath, mode='w')

        # Create level zero groups
        chromosomes_group = out_f.create_group("chromosomes")

        # Prepare to fill in chroms dataset
        chromosomes = in_bins_df["chr_name"].unique().tolist()
        num_chromosomes = len(chromosomes)

        in_chrom_ends_df = in_bins_df.drop_duplicates(subset=['chr_name'], keep='last')
        in_chrom_ends_df = in_chrom_ends_df.set_index("chr_name")

        chroms_length_arr = np.array([ in_chrom_ends_df.at[x, "chr_end"] for x in chromosomes ], dtype="i8")
        chroms_cumsum_arr = np.concatenate((np.array([0]), np.cumsum(chroms_length_arr)))

        chromosomes_set = set(chromosomes)
        chrom_name_to_length = dict(zip(chromosomes, chroms_length_arr))
        chrom_name_to_cumsum = dict(zip(chromosomes, chroms_cumsum_arr))

        
        # Prepare to fill in resolutions dataset
        resolutions = [ starting_resolution*(2**x) for x in range(16)]
        
        # Create each chromosome dataset.
        for chr_name, chr_len in chrom_name_to_length.items():
            chr_group = chromosomes_group.create_group(chr_name)
            # Create each resolution group.
            for resolution in resolutions:
                chr_shape = (num_samples, math.ceil(chr_len / resolution))
                chr_group.create_dataset(str(resolution), shape=chr_shape, dtype="f4", fill_value=np.nan, compressor=compressor)
        
        # Fill in data for each cluster.
        in_clusters_df["cluster"] = in_clusters_df["cluster"].astype(str)
        cluster_ids = in_clusters_df["cluster"].unique().tolist()
        cluster_ids.sort(key=int)

        cluster_profiles = {}
        for cluster_index, cluster_id in enumerate(cluster_ids):
            cluster_df = in_clusters_df.loc[in_clusters_df["cluster"] == cluster_id]
            cluster_cell_ids = cluster_df.index.values.tolist()
            cluster_num_cells = len(chrom_num_bins)
            cluster_cell_indices = in_barcodes_df.loc[in_barcodes_df[0].isin(cluster_cell_ids)].index.values


            for chrom_name in chromosomes:
                chrom_len = chrom_name_to_length[chrom_name]
                chrom_num_bins = math.ceil(chrom_len / bin_size)
                chrom_indices = in_bins_df.loc[in_bins_df["chr_name"] == chrom_name].index.values
                
                cluster_cell_by_bin_mtx = in_mtx[cluster_cell_indices, chrom_indices]

                cluster_profiles[cluster_id] = cluster_cell_by_bin_mtx


            # Fill in data for each resolution of a bigwig file.
            for resolution in resolutions:
                # Fill in data for each chromosome of a resolution of a bigwig file.
                for chr_name in matching_chromosomes:
                    chr_len = chrom_name_to_length[chr_name]
                    chr_shape = (num_samples, math.ceil(chr_len / resolution))
                    # Group every `resolution` values together and take sum.
                    arr = np.reshape(values, (-1, resolution)).sum(axis=-1)
                    chromosomes_group[chr_name][str(resolution)][cluster_index,:] = arr
        
        max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print(max_mem)
        
        # Append metadata to the top resolution row_infos attribute.
        row_infos = []
        for metadata_index, metadata_file in enumerate(input_metadata_files):
            with open(metadata_file) as mf:
                metadata_json = json.load(mf)
            row_info = metadata_json_to_row_info(metadata_json)
            row_infos.append(row_info)
        
        # f.attrs should contain all tileset_info properties
        # For zarr, more attributes are used here to allow "serverless"
        f.attrs['row_infos'] = row_infos
        f.attrs['resolutions'] = sorted(resolutions, reverse=True)
        f.attrs['shape'] = [ num_samples, 256 ]
        f.attrs['name'] = name
        f.attrs['coordSystem'] = name_to_coordsystem(name)
        
        # https://github.com/zarr-developers/zarr-specs/issues/50
        f.attrs['multiscales'] = [
            {
                "version": "0.1",
                "name": chr_name,
                "datasets": [
                    { "path": f"chromosomes/{chr_name}/{resolution}" }
                    for resolution in sorted(resolutions, reverse=True)
                ],
                "type": "zarr-multivec",
                "metadata": {
                    "chromoffset": int(chrom_name_to_cumsum[chr_name]),
                    "chromsize": int(chr_len),
                }
            }
            for (chr_name, chr_len) in list(zip(chromosomes, chroms_length_arr))
        ]


        return

    def get_genomic_profiles(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
        
        zarr_tempdir = self.tempdir
        zarr_filepath = join(zarr_tempdir, 'profiles.zarr')

        self._create_genomic_multivec_zarr(zarr_filepath)

        if zarr_tempdir is not None:
            obj_routes = [
                Mount(self._get_route(dataset_uid, obj_i, "genomic"),
                    app=StaticFiles(directory=os.path.dirname(zarr_filepath), html=False, check_dir=False)),
            ]

            obj_file_defs = [
                {
                    "type": dt.GENOMIC_PROFILES.value,
                    "fileType": ft.GENOMIC_PROFILES_ZARR.value,
                    "url": self._get_url(port, dataset_uid, obj_i, "genomic/profiles.zarr")
                }
            ]

        return obj_file_defs, obj_routes