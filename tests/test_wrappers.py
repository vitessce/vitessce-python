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
    JsonWrapper,
    MultivecZarrWrapper,
    ImageOmeTiffWrapper,
    ObsSegmentationsOmeTiffWrapper,
    ImageOmeZarrWrapper,
    ObsSegmentationsOmeZarrWrapper,
)

from vitessce.wrappers import SpatialDataWrapper, file_path_to_url_path


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

    def test_image_ome_tiff(self):
        w = ImageOmeTiffWrapper(img_path=data_path / 'test.ome.tif')
        w.local_img_uid = 'test.ome.tif'
        w.local_offsets_uid = 'test.offsets.json'

        file_def_creator = w.make_raster_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'image.ome-tiff',
            'url': 'http://localhost:8000/A/0/test.ome.tif',
            'options': {
                'offsetsUrl': 'http://localhost:8000/A/0/test.offsets.json'
            }
        })

    def test_image_ome_zarr(self):
        w = ImageOmeZarrWrapper(img_path=data_path / 'test.ome.zarr')
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

    def test_obs_segmentations_ome_tiff(self):
        w = ObsSegmentationsOmeTiffWrapper(img_path=data_path / 'test.ome.tif')
        w.local_img_uid = 'test.ome.tif'
        w.local_offsets_uid = 'test.offsets.json'

        file_def_creator = w.make_raster_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsSegmentations.ome-tiff',
            'url': 'http://localhost:8000/A/0/test.ome.tif',
            'options': {
                'offsetsUrl': 'http://localhost:8000/A/0/test.offsets.json'
            }
        })

    def test_obs_segmentations_ome_zarr(self):
        w = ObsSegmentationsOmeZarrWrapper(img_path=data_path / 'test.ome.zarr')
        w.local_dir_uid = 'test.ome.zarr'

        file_def_creator = w.make_image_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsSegmentations.ome-zarr',
            'url': 'http://localhost:8000/A/0/test.ome.zarr'
        })

    def test_anndata(self):
        adata_path = data_path / 'test.h5ad.zarr'
        w = AnnDataWrapper(adata_path,
                           obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'],
                           obs_labels_names=['Cell Label'], obs_labels_paths=['obs/CellLabel'],
                           obs_embedding_paths=['obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.local_dir_uid = 'anndata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.zarr', 'url': 'http://localhost:8000/A/0/anndata.zarr',
                                    'options': {
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'embeddingType': 'UMAP', 'dims': [0, 1]}],
                                        'obsSets': [{'path': 'obs/CellType', 'name': 'Cell Type'}],
                                        'obsLabels': [{'path': 'obs/CellLabel', 'obsLabelsType': 'Cell Label'}]
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
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'dims': [0, 1], 'embeddingType': 'UMAP'}],
                                        'obsSets': [{'name': 'Cell Type', 'path': 'obs/CellType'}]
                                    }})

    def test_anndata_with_base_dir_no_names(self):
        adata_path = 'test.h5ad.zarr'
        w = AnnDataWrapper(adata_path, obs_set_paths=['obs/CellType'], obs_embedding_paths=[
                           'obsm/X_umap'])
        w.base_dir = data_path
        w.local_dir_uid = 'anndata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.zarr', 'url': 'http://localhost:8000/test.h5ad.zarr',
                                    'options': {
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'dims': [0, 1], 'embeddingType': 'X_umap'}],
                                        'obsSets': [{'name': 'CellType', 'path': 'obs/CellType'}]
                                    }})

    def test_anndata_with_h5ad_and_ref_json(self):
        adata_path = data_path / 'test.h5ad'
        ref_json_path = data_path / 'test.h5ad.ref.json'
        w = AnnDataWrapper(adata_path, ref_path=ref_json_path,
                           obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'],
                           obs_labels_names=['Cell Label'], obs_labels_paths=['obs/CellLabel'],
                           obs_embedding_paths=['obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.local_file_uid = 'anndata.h5ad'
        w.local_ref_uid = 'anndata.reference.json'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.h5ad', 'url': 'http://localhost:8000/A/0/anndata.h5ad',
                                    'options': {
                                        'refSpecUrl': 'http://localhost:8000/A/0/anndata.reference.json',
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'embeddingType': 'UMAP', 'dims': [0, 1]}],
                                        'obsSets': [{'path': 'obs/CellType', 'name': 'Cell Type'}],
                                        'obsLabels': [{'path': 'obs/CellLabel', 'obsLabelsType': 'Cell Label'}]
                                    }})

    def test_anndata_with_h5ad_and_ref_json_with_base_dir(self):
        adata_path = 'test.h5ad'
        ref_json_path = 'test.h5ad.ref.json'
        w = AnnDataWrapper(adata_path, ref_path=ref_json_path,
                           obs_set_paths=['obs/CellType'], obs_set_names=['Cell Type'],
                           obs_labels_names=['Cell Label'], obs_labels_paths=['obs/CellLabel'],
                           obs_embedding_paths=['obsm/X_umap'], obs_embedding_names=['UMAP'])
        w.base_dir = data_path
        w.local_file_uid = 'anndata.h5ad'
        w.local_ref_uid = 'anndata.reference.json'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {'fileType': 'anndata.h5ad', 'url': 'http://localhost:8000/test.h5ad',
                                    'options': {
                                        'refSpecUrl': 'http://localhost:8000/test.h5ad.ref.json',
                                        'obsEmbedding': [{'path': 'obsm/X_umap', 'embeddingType': 'UMAP', 'dims': [0, 1]}],
                                        'obsSets': [{'path': 'obs/CellType', 'name': 'Cell Type'}],
                                        'obsLabels': [{'path': 'obs/CellLabel', 'obsLabelsType': 'Cell Label'}]
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

    def test_json(self):
        w = JsonWrapper(
            json_path=data_path / 'test.segmentations.json',
            data_type="obsSegmentations",
            coordination_values={
                "obsType": "nucleus"
            }
        )
        w.local_json_uid = 'test_uid'

        file_def_creator = w.make_json_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsSegmentations.json',
            'url': 'http://localhost:8000/A/0/test_uid',
            'coordinationValues': {
                "obsType": "nucleus"
            }
        })

    def test_json_with_base_dir(self):
        w = JsonWrapper(
            json_path='test.segmentations.json',
            data_type="obsSegmentations",
            coordination_values={
                "obsType": "nucleus"
            }
        )
        w.base_dir = data_path
        w.local_json_uid = 'test_uid'

        file_def_creator = w.make_json_file_def_creator(
            "A",
            "0"
        )
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'obsSegmentations.json',
            'url': 'http://localhost:8000/test.segmentations.json',
            'coordinationValues': {
                "obsType": "nucleus"
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

    def test_spatial_data_with_base_dir(self):

        spatial_data_path = 'test.spatialdata.zarr'
        w = SpatialDataWrapper(
            sdata_path=spatial_data_path,
            image_path="images/picture",
            obs_set_paths=['obs/CellType'],
            obs_set_names=['Cell Type'],
            obs_embedding_paths=['obsm/X_umap'],
            obs_embedding_names=['UMAP']
        )
        w.base_dir = data_path
        w.local_dir_uid = 'spatialdata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertEqual(file_def, {
            'fileType': 'spatialdata.zarr',
            'url': 'http://localhost:8000/test.spatialdata.zarr',
            'options': {
                'obsSets': {
                    'obsSets': [{'name': 'Cell Type', 'path': 'obs/CellType'}],
                    'tablePath': 'tables/table'
                },
                'image': {'path': 'images/picture'}
            }})

    def test_spatial_data_with_base_dir_2(self):
        spatial_data_path = 'test.spatialdata.zarr'
        w = SpatialDataWrapper(
            sdata_path=spatial_data_path,
            image_path='images/CytAssist_FFPE_Human_Breast_Cancer_full_image',
            coordinate_system='aligned',
            region='CytAssist_FFPE_Human_Breast_Cancer',
            obs_feature_matrix_path='tables/table/X',
            obs_spots_path='shapes/CytAssist_FFPE_Human_Breast_Cancer',
            table_path='tables/table',
            coordination_values={
                "obsType": "spot"
            }
        )
        w.base_dir = data_path
        w.local_dir_uid = 'spatialdata.zarr'

        file_def_creator = w.make_file_def_creator('A', 0)
        file_def = file_def_creator('http://localhost:8000')
        self.assertDictEqual(file_def, {
            'fileType': 'spatialdata.zarr',
            'url': 'http://localhost:8000/test.spatialdata.zarr',
            'options': {
                'image': {
                    'path': 'images/CytAssist_FFPE_Human_Breast_Cancer_full_image',
                    'coordinateSystem': 'aligned',
                },
                'obsFeatureMatrix': {
                    'path': 'tables/table/X',
                    'region': 'CytAssist_FFPE_Human_Breast_Cancer'
                },
                'obsSpots': {
                    'path': 'shapes/CytAssist_FFPE_Human_Breast_Cancer',
                    'tablePath': 'tables/table',
                    'region': 'CytAssist_FFPE_Human_Breast_Cancer',
                    'coordinateSystem': 'aligned',
                }
            },
            'coordinationValues': {
                "obsType": "spot"
            }
        })
