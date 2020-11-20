import os
import tempfile

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
        offsets_dir_path, offsets_url = (None, None) if self.offsets_path is None else (self._get_offsets_dir(), self._get_url(port, dataset_uid, obj_i, os.path.join("raster_offsets", self._get_offsets_filename())))

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
        try:
            import zarr
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
        except ImportError:
            pass

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
    
    def _create_exp_matrix_zarr(self):
        
        try:
            import zarr
            from numcodecs import Zlib
            from scipy.sparse import csr_matrix

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

            zarr_tempdir = self.tempdir
            zarr_filepath = os.path.join(zarr_tempdir, 'matrix.zarr')
        
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
        except ImportError:
            zarr_tempdir = None
            zarr_filepath = None
        
        return zarr_tempdir, zarr_filepath


    def get_cells(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
        try:
            import anndata
            if type(self.adata) == anndata.AnnData:
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
        except ImportError:
            pass

        return obj_file_defs, obj_routes

    def get_cell_sets(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []
        try:
            import anndata
            if type(self.adata) == anndata.AnnData:
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
        except ImportError:
            pass

        return obj_file_defs, obj_routes
    
    def get_expression_matrix(self, port, dataset_uid, obj_i):
        obj_routes = []
        obj_file_defs = []

        zarr_tempdir, zarr_filepath = self._create_exp_matrix_zarr()

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
