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

        raster_file_def_creator = w.make_raster_file_def_creator(
            "A",
            "0"
        )
        raster_json = raster_file_def_creator('http://localhost:8000')
        self.assertEqual(raster_json, {
            'type': 'raster',
            'fileType': 'raster.json',
            'options': {
                'schemaVersion': '0.0.2',
                'images': [
                    {
                        'name': 'Test',
                        'type': 'ome-tiff',
                        'url': 'http://localhost:8000/A/0/test.ome.tif',
                        'metadata': {
                            'omeTiffOffsetsUrl': 'http://localhost:8000/A/0/test.offsets.json'
                        }
                    }
                ],
            }
        })
    
    def test_anndata(self):
        adata = read_h5ad(join('data', 'test.h5ad'))
        w = AnnDataWrapper(adata, cell_set_obs=['CellType'])

        cells_creator = w.make_cells_file_def_creator('A', 0)
        cells = cells_creator( 'http://localhost:8000')
        self.assertEqual(cells, {'type': 'cells', 'fileType': 'anndata-cells.zarr', 'url': 'http://localhost:8000/A/0/anndata.zarr', 'options': {'factors': ['obs/CellType']}})

        cell_sets_creator = w.make_cell_sets_file_def_creator('A', 0)
        cell_sets = cell_sets_creator( 'http://localhost:8000')
        self.assertEqual(cell_sets, {'type': 'cell-sets', 'fileType': 'anndata-cell-sets.zarr', 'url': 'http://localhost:8000/A/0/anndata.zarr', 'options': [{'groupName': 'CellType', 'setName': 'obs/CellType'}]})

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
        