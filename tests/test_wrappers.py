import unittest
from os.path import join

import zarr
from anndata import read_h5ad
from scipy.io import mmread
import pandas as pd
import numpy as np

from create_test_data import (
    create_test_anndata_file,
    create_test_loom_file,
    create_test_ometiff_file,
    create_test_omezarr_store,
    create_test_snaptools_files,
)

from vitessce import (
    OmeTiffWrapper,
    OmeZarrWrapper,
    AnnDataWrapper,
    SnapWrapper,
)

class TestWrappers(unittest.TestCase):

    def setUp(self):
        create_test_anndata_file(join('data', 'test.h5ad'))
        create_test_loom_file(join('data', 'test.loom'))
        create_test_ometiff_file(join('data', 'test.ome.tif'))
        create_test_omezarr_store(join('data', 'test.ome.zarr'))
        create_test_snaptools_files(
            join('data', 'test.snap.mtx'),
            join('data', 'test.snap.bins.txt'),
            join('data', 'test.snap.barcodes.txt'),
            join('data', 'test.snap.clusters.csv'),
        )

    def test_ome_tiff(self):
        w = OmeTiffWrapper(img_path="data/test.ome.tif", name="Test")

        raster_json = w.create_raster_json(
            "http://localhost:8000/A/0/test.ome.tif",
            "http://localhost:8000/A/0/test.offsets.json"
        )

        self.assertEqual(raster_json, {
            'images': [
                {
                    'metadata': {
                        'omeTiffOffsetsUrl': 'http://localhost:8000/A/0/test.offsets.json'
                    },
                    'name': 'Test',
                    'type': 'ome-tiff',
                    'url': 'http://localhost:8000/A/0/test.ome.tif'
                }
            ],
            'schemaVersion': '0.0.2'
        })

        obj_file_defs, obj_routes = w.get_raster("http://localhost:8000", 'A', 0)

        self.assertEqual(obj_file_defs, [
            {
                'type': 'raster',
                'fileType': 'raster.json',
                'options': raster_json
            }
        ])
    
    def test_omezarr_store(self):
        z = zarr.open('data/test.ome.zarr')
        w = OmeZarrWrapper(z)

        raster_json = w.create_raster_json(
            "http://localhost:8000/raster_img"
        )
        
        # TODO
        # self.assertEqual(raster_json, {})

        obj_file_defs, obj_routes = w.get_raster("http://localhost:8000", 'A', 0)
        self.assertEqual(obj_file_defs, [
            {
                'fileType': 'raster.json',
                'type': 'raster',
                'url': 'http://localhost:8000/A/0/raster'
            }
        ])
    
    def test_base_url(self):
        z = zarr.open('data/test.ome.zarr')
        w = OmeZarrWrapper(z)

        raster_json = w.create_raster_json(
            "https://example.com/raster_img"
        )
        
        # TODO
        # self.assertEqual(raster_json, {})

        obj_file_defs, obj_routes = w.get_raster("https://example.com", 'A', 0)
        self.assertEqual(obj_file_defs, [
            {
                'fileType': 'raster.json',
                'type': 'raster',
                'url': 'https://example.com/A/0/raster'
            }
        ])
    
    def test_anndata(self):
        adata = read_h5ad(join('data', 'test.h5ad'))
        w = AnnDataWrapper(adata, cell_set_obs_cols=['CellType'])

        cells_json = w.create_cells_json()
        cell_sets_json = w.create_cell_sets_json()

        obj_file_defs, obj_routes = w.get_cells("http://localhost:8000", 'A', 0)
        self.assertEqual(obj_file_defs, [{'type': 'cells', 'fileType': 'cells.json', 'url': 'http://localhost:8000/A/0/cells'}])

        obj_file_defs, obj_routes = w.get_cell_sets("http://localhost:8000", 'A', 0)
        self.assertEqual(obj_file_defs, [{'type': 'cell-sets', 'fileType': 'cell-sets.json', 'url': 'http://localhost:8000/A/0/cell-sets'}])

    def test_snaptools(self):
        mtx = mmread(join('data', 'test.snap.mtx'))
        barcodes_df = pd.read_csv(join('data', 'test.snap.barcodes.txt'), header=None)
        bins_df = pd.read_csv(join('data', 'test.snap.bins.txt'), header=None)
        clusters_df = pd.read_csv(join('data', 'test.snap.clusters.csv'), index_col=0)

        zarr_filepath = join('data', 'test_out.snap.multivec.zarr')

        w = SnapWrapper(mtx, barcodes_df, bins_df, clusters_df)
        w.create_genomic_multivec_zarr(zarr_filepath)

        z = zarr.open(zarr_filepath, mode='r')

        self.assertEqual(z['chromosomes/chr1/5000'].shape, (4, 49792))
        self.assertEqual(z['chromosomes/chr1/5000'][:,0].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][:,1].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][:,2].sum(), 7)
        self.assertEqual(z['chromosomes/chr1/5000'][:,3].sum(), 7)
        self.assertEqual(z['chromosomes/chr1/5000'][:,4].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][0,:].sum(), 17)
        self.assertEqual(z['chromosomes/chr1/10000'][0,:].sum(), 17)
        self.assertEqual(z['chromosomes/chr1/5000'][:,2].sum(), 7)
        self.assertEqual(z['chromosomes/chr2/5000'][:,0].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/5000'][:,1].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/10000'][:,0].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/5000'][:,2].sum(), 4)
        self.assertEqual(z['chromosomes/chr2/5000'][:,3].sum(), 9)
        self.assertEqual(z['chromosomes/chr2/10000'][:,1].sum(), 13)
        self.assertEqual(z['chromosomes/chr3/5000'][:,3].sum(), 9)
        self.assertEqual(z['chromosomes/chr3/5000'][:].sum(), 9)
        self.assertEqual(z['chromosomes/chr18/5000'][:].sum(), 8)

        cells_json = w.create_cells_json()
        self.assertEqual(len(cells_json), 6)
        self.assertEqual(cells_json['AAACATCGAGTACAAGACAGCAGA'], { 'mappings': { 'UMAP': [4.43, 1.64] } })
        