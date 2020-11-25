import unittest
from os.path import join

import zarr
from anndata import read_h5ad
from scipy.io import mmread
import pandas as pd
import numpy as np

from vitessce import (
    OmeTiffWrapper,
    ZarrDirectoryStoreWrapper,
    AnnDataWrapper,
    LoomWrapper,
    SnapToolsWrapper,
)

class TestWrappers(unittest.TestCase):

    def test_ome_tiff(self):
        w = OmeTiffWrapper("data/test.ome.tif", offsets_path="data/offsets.json", name="Test")

        raster_json = w._create_raster_json(
            "http://localhost:8000/raster_img",
            "http://localhost:8000/raster_offsets/offsets.json"
        )

        self.assertEqual(raster_json, {
            'images': [
                {
                    'metadata': {
                        'omeTiffOffsetsUrl': 'http://localhost:8000/raster_offsets/offsets.json'
                    },
                    'name': 'Test',
                    'type': 'ome-tiff',
                    'url': 'http://localhost:8000/raster_img'
                }
            ],
            'schemaVersion': '0.0.2'
        })

        obj_file_defs, obj_routes = w.get_raster(8000, 'A', 0)

        self.assertEqual(obj_file_defs, [
            {
                'fileType': 'raster.json',
                'type': 'raster',
                'url': 'http://localhost:8000/A/0/raster'
            }
        ])
    
    def test_zarr_dir_store(self):
        z = zarr.open('data/test.ome.zarr')
        w = ZarrDirectoryStoreWrapper(z)

        raster_json = w._create_raster_json(
            "http://localhost:8000/raster_img"
        )
        
        # TODO
        # self.assertEqual(raster_json, {})

        obj_file_defs, obj_routes = w.get_raster(8000, 'A', 0)
        self.assertEqual(obj_file_defs, [
            {
                'fileType': 'raster.json',
                'type': 'raster',
                'url': 'http://localhost:8000/A/0/raster'
            }
        ])
    
    def test_anndata(self):
        adata = read_h5ad(join('data', 'habib17.processed.h5ad'))
        w = AnnDataWrapper(adata)

        cells_json = w._create_cells_json()
        cell_sets_json = w._create_cell_sets_json()

        obj_file_defs, obj_routes = w.get_cells(8000, 'A', 0)
        self.assertEqual(obj_file_defs, [{'type': 'cells', 'fileType': 'cells.json', 'url': 'http://localhost:8000/A/0/cells'}])

        obj_file_defs, obj_routes = w.get_cell_sets(8000, 'A', 0)
        self.assertEqual(obj_file_defs, [{'type': 'cell-sets', 'fileType': 'cell-sets.json', 'url': 'http://localhost:8000/A/0/cell-sets'}])

    def test_snaptools(self):
        mtx = mmread(join('data', 'filtered_cell_by_bin.mtx'))
        barcodes_df = pd.read_csv(join('data', 'barcodes.txt'), header=None)
        bins_df = pd.read_csv(join('data', 'bins.txt'), header=None)
        clusters_df = pd.read_csv(join('data', 'umap_coords_clusters.csv'), index_col=0)

        print(bins_df.head())

        zarr_filepath = join('data', 'test_snaptools.zarr')

        w = SnapToolsWrapper(mtx, barcodes_df, bins_df, clusters_df)
        w._create_genomic_multivec_zarr(zarr_filepath)

        z = zarr.open(zarr_filepath, mode='r')
        print(z)
        