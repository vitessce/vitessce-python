import unittest

from anndata import read_h5ad
from pathlib import Path

from .create_test_data import (
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
    CsvWrapper,
    MultivecZarrWrapper,
)


data_path = Path('tests/data')


class TestWrappers(unittest.TestCase):

    def setUp(self):
        create_test_anndata_file(data_path / 'test.h5ad')
        create_test_loom_file(data_path / 'test.loom')
        create_test_ometiff_file(data_path / 'test.ome.tif')
        create_test_omezarr_store(data_path / 'test.ome.zarr')
        create_test_snaptools_files(
            data_path / 'test.snap.mtx',
            data_path / 'test.snap.bins.txt',
            data_path / 'test.snap.barcodes.txt',
            data_path / 'test.snap.clusters.csv',
        )

    def test_ome_tiff(self):
        w = OmeTiffWrapper(img_path=data_path / 'test.ome.tif', name="Test")
        w.local_img_uid = 'test.ome.tif'
        w.local_offsets_uid = 'test.offsets.json'

        file_def_creator = w.make_raster_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'raster.json',
            'options': {
                'schemaVersion': '0.0.2',
                'images': [
                    {
                        'name': 'Test',
                        'type': 'ome-tiff',
                        'url': 'http://localhost:8000/A/0/test.ome.tif',
                        'metadata': {
                            'isBitmask': False,
                            'omeTiffOffsetsUrl': 'http://localhost:8000/A/0/test.offsets.json'
                        }
                    }
                ],
            }
        })

    def test_ome_zarr(self):
        w = OmeZarrWrapper(img_path=data_path / 'test.ome.zarr')
        w.local_dir_uid = 'test.ome.zarr'

        file_def_creator = w.make_image_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'image.ome-zarr',
            'url': 'http://localhost:8000/A/0/test.ome.zarr'
        })

    def test_anndata(self):
        adata = read_h5ad(data_path / 'test.h5ad')
        w = AnnDataWrapper(adata, obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'], obs_embedding_paths=[
                           'obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.local_dir_uid = 'anndata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.zarr', 'url': 'http://localhost:8000/A/0/anndata.zarr',
                                    'options': {
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'embeddingType': 'UMAP', 'dims': [0, 1]}],
                                        'obsSets': [{'path': 'obs/CellType', 'name': 'Cell Type'}]
                                    }})

    def test_csv(self):
        w = CsvWrapper(
            csv_path=data_path / 'test.umap.csv',
            data_type="obsEmbedding",
            options={
                "obsIndex": "index",
                "obsEmbedding": ["UMAP_1", "UMAP_2"]
            },
            coordination_values={
                "embeddingType": "UMAP"
            }
        )
        w.local_csv_uid = 'test_uid'

        file_def_creator = w.make_csv_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsEmbedding.csv',
            'url': 'http://localhost:8000/A/0/test_uid',
            'options': {
                "obsIndex": "index",
                "obsEmbedding": ["UMAP_1", "UMAP_2"]
            },
            'coordinationValues': {
                "embeddingType": "UMAP"
            }
        })

    def test_multivec_zarr(self):
        zarr_filepath = data_path / 'test_out.snap.multivec.zarr'

        w = MultivecZarrWrapper(zarr_path=zarr_filepath)
        w.local_dir_uid = 'test_uid'

        file_def_creator = w.make_genomic_profiles_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'genomic-profiles.zarr',
            'url': 'http://localhost:8000/A/0/test_uid',
        })
