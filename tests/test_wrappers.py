import unittest

from pathlib import Path, PurePosixPath, PureWindowsPath

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

from vitessce.wrappers import file_path_to_url_path


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

    def test_file_path_to_url_path(self):
        posix_with_slash = file_path_to_url_path("tests/data/test.snap.mtx", path_class=PurePosixPath)
        self.assertEqual(posix_with_slash, "/tests/data/test.snap.mtx")

        posix_without_slash = file_path_to_url_path("tests/data/test.snap.mtx", prepend_slash=False, path_class=PurePosixPath)
        self.assertEqual(posix_without_slash, "tests/data/test.snap.mtx")

        posix_with_dot_and_slash = file_path_to_url_path("./tests/data/test.snap.mtx", path_class=PurePosixPath)
        self.assertEqual(posix_with_dot_and_slash, "/tests/data/test.snap.mtx")

        posix_with_dot_without_slash = file_path_to_url_path("./tests/data/test.snap.mtx", prepend_slash=False, path_class=PurePosixPath)
        self.assertEqual(posix_with_dot_without_slash, "tests/data/test.snap.mtx")

        windows_with_slash = file_path_to_url_path("tests\\data\\test.snap.mtx", path_class=PureWindowsPath)
        self.assertEqual(windows_with_slash, "/tests/data/test.snap.mtx")

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

    def test_ome_tiff_with_base_dir(self):
        w = OmeTiffWrapper(img_path='test.ome.tif', name="Test")
        w.base_dir = str(data_path)
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
                        'url': 'http://localhost:8000/test.ome.tif',
                        'metadata': {
                            'isBitmask': False,
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

    def test_ome_zarr_with_base_dir(self):
        w = OmeZarrWrapper(img_path='test.ome.zarr')
        w.base_dir = str(data_path)
        w.local_dir_uid = 'test.ome.zarr'

        file_def_creator = w.make_image_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'image.ome-zarr',
            'url': 'http://localhost:8000/test.ome.zarr'
        })

    def test_anndata(self):
        adata_path = data_path / 'test.h5ad.zarr'
        w = AnnDataWrapper(adata_path, obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'], obs_embedding_paths=[
                           'obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.local_dir_uid = 'anndata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.zarr', 'url': 'http://localhost:8000/A/0/anndata.zarr',
                                    'options': {
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'embeddingType': 'UMAP', 'dims': [0, 1]}],
                                        'obsSets': [{'path': 'obs/CellType', 'name': 'Cell Type'}]
                                    }})

    def test_anndata_with_base_dir(self):
        adata_path = 'test.h5ad.zarr'
        w = AnnDataWrapper(adata_path, obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'], obs_embedding_paths=[
                           'obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.base_dir = data_path
        w.local_dir_uid = 'anndata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.zarr', 'url': 'http://localhost:8000/test.h5ad.zarr',
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

    def test_csv_with_base_dir(self):
        w = CsvWrapper(
            csv_path='test.umap.csv',
            data_type="obsEmbedding",
            options={
                "obsIndex": "index",
                "obsEmbedding": ["UMAP_1", "UMAP_2"]
            },
            coordination_values={
                "embeddingType": "UMAP"
            }
        )
        w.base_dir = data_path
        w.local_csv_uid = 'test_uid'

        file_def_creator = w.make_csv_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsEmbedding.csv',
            'url': 'http://localhost:8000/test.umap.csv',
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

    def test_multivec_zarr_with_base_dir(self):
        zarr_filepath = 'test_out.snap.multivec.zarr'

        w = MultivecZarrWrapper(zarr_path=zarr_filepath)
        w.base_dir = data_path
        w.local_dir_uid = 'test_uid'

        file_def_creator = w.make_genomic_profiles_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'genomic-profiles.zarr',
            'url': 'http://localhost:8000/test_out.snap.multivec.zarr',
        })
