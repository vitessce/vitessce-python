import pytest
import ast

from vitessce import (  # noqa: F401
    VitessceConfig,
    CoordinationType as ct,
    ViewType as cm,
    DataType as dt,
    FileType as ft,
    hconcat,
    vconcat,
    AbstractWrapper,
    make_repr,
    CoordinationLevel as CL,
    AnnDataWrapper,

    # Neither of these is in the source code, but they do appear in code which is eval'd.
    VitessceChainableConfig,
    VitessceConfigDatasetFile,
)


class MockArtifactPath:
    def __init__(self, url):
        self.url = url

    def to_url(self):
        return self.url


class MockArtifact:
    def __init__(self, name, url):
        self.name = name
        self.path = MockArtifactPath(url)


def test_config_creation():
    vc = VitessceConfig(schema_version="1.0.15")
    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [],
        "coordinationSpace": {},
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_add_dataset():
    vc = VitessceConfig(schema_version="1.0.15")
    vc.add_dataset(name='My Dataset')

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_add_anndata_url():
    vc = VitessceConfig(schema_version="1.0.15")
    vc.add_dataset(name='My Dataset').add_object(
        AnnDataWrapper(
            adata_url="http://example.com/adata.h5ad.zarr",
            obs_set_paths=["obs/louvain"],
        )
    )

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': [
                    {
                        "fileType": "anndata.zarr",
                        "url": "http://example.com/adata.h5ad.zarr",
                        "options": {
                            "obsSets": [
                                {
                                    "name": "louvain",
                                    "path": "obs/louvain",
                                }
                            ]
                        }
                    }
                ]
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_add_anndata_artifact():
    vc = VitessceConfig(schema_version="1.0.15")
    vc.add_dataset(name='My Dataset').add_object(
        AnnDataWrapper(
            adata_artifact=MockArtifact("My anndata artifact", "http://example.com/adata.h5ad.zarr"),
            obs_set_paths=["obs/louvain"],
        )
    )

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': [
                    {
                        "fileType": "anndata.zarr",
                        "url": "http://example.com/adata.h5ad.zarr",
                        "options": {
                            "obsSets": [
                                {
                                    "name": "louvain",
                                    "path": "obs/louvain",
                                }
                            ]
                        }
                    }
                ]
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [],
        "initStrategy": "auto"
    }

    vc_artifacts = vc.get_artifacts()
    assert list(vc_artifacts.keys()) == ["http://example.com/adata.h5ad.zarr"]
    assert vc_artifacts["http://example.com/adata.h5ad.zarr"].name == "My anndata artifact"


def test_config_add_dataset_add_files():
    vc = VitessceConfig(schema_version="1.0.15")
    vc.add_dataset(name='My Chained Dataset').add_file(
        url="http://example.com/cells.json",
        file_type=ft.CELLS_JSON,
        coordination_values={"obsType": "cell"},
    ).add_file(
        url="http://example.com/cell_sets.json",
        file_type=ft.CELL_SETS_JSON,
        coordination_values={"obsType": "cell"},
    )

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Chained Dataset',
                'files': [
                    {
                        'url': 'http://example.com/cells.json',
                        'fileType': 'cells.json',
                        'coordinationValues': {
                            'obsType': 'cell'
                        }
                    },
                    {
                        'url': 'http://example.com/cell_sets.json',
                        'fileType': 'cell-sets.json',
                        'coordinationValues': {
                            'obsType': 'cell'
                        }
                    }
                ]
            },
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_add_spatial_view():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')

    vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'h': 1,
                'w': 1,
                'x': 0,
                'y': 0
            }
        ],
        "initStrategy": "auto"
    }


def test_config_add_scatterplot_view_with_mapping():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')

    vc.add_view(
        cm.SCATTERPLOT, dataset=my_dataset, mapping="X_umap")

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
            'embeddingType': {
                'A': 'X_umap'
            }
        },
        "layout": [
            {
                'component': 'scatterplot',
                'coordinationScopes': {
                    'dataset': 'A',
                    'embeddingType': 'A'
                },
                'h': 1,
                'w': 1,
                'x': 0,
                'y': 0
            }
        ],
        "initStrategy": "auto"
    }


def test_config_add_scatterplot_view_with_embedding_coordinations():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')

    my_view = vc.add_view(cm.SCATTERPLOT, dataset=my_dataset)

    et_scope, ez_scope, ex_scope, ey_scope = vc.add_coordination(
        ct.EMBEDDING_TYPE, ct.EMBEDDING_ZOOM, ct.EMBEDDING_TARGET_X, ct.EMBEDDING_TARGET_Y)
    my_view.use_coordination(et_scope, ez_scope, ex_scope, ey_scope)

    et_scope.set_value("X_pca")
    ez_scope.set_value(2)
    ex_scope.set_value(10)
    ey_scope.set_value(11)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
            'embeddingType': {
                'A': 'X_pca'
            },
            'embeddingZoom': {
                'A': 2
            },
            'embeddingTargetX': {
                'A': 10
            },
            'embeddingTargetY': {
                'A': 11
            },
        },
        "layout": [
            {
                'component': 'scatterplot',
                'coordinationScopes': {
                    'dataset': 'A',
                    'embeddingType': 'A',
                    'embeddingZoom': 'A',
                    'embeddingTargetX': 'A',
                    'embeddingTargetY': 'A',
                },
                'h': 1,
                'w': 1,
                'x': 0,
                'y': 0
            }
        ],
        "initStrategy": "auto"
    }


def test_config_add_dataset_add_objects():
    vc = VitessceConfig(schema_version="1.0.15")

    class MockWrapperA(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i, base_dir=None):
            def get_molecules(base_url):
                return {
                    "url": f"{base_url}/molecules",
                    "fileType": "molecules.json",
                    "coordinationValues": {
                        "obsType": "molecule"
                    }
                }

            def get_cells(base_url):
                return {
                    "url": f"{base_url}/cells",
                    "fileType": "cells.json"
                }
            self.file_def_creators += [get_molecules, get_cells]

    class MockWrapperB(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i, base_dir=None):
            def get_cell_sets(base_url):
                return {
                    "url": f"{base_url}/cell-sets",
                    "fileType": "cell-sets.json"
                }
            self.file_def_creators += [get_cell_sets]

    vc.add_dataset(name='My Object Dataset').add_object(
        obj=MockWrapperA("Experiment A")
    ).add_object(
        obj=MockWrapperB("Experiment B")
    )

    vc_dict = vc.to_dict(base_url="http://localhost:8000")

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Object Dataset',
                'files': [
                    {
                        "url": "http://localhost:8000/molecules",
                        "fileType": "molecules.json",
                        "coordinationValues": {
                            "obsType": "molecule"
                        }
                    },
                    {
                        "url": "http://localhost:8000/cells",
                        "fileType": "cells.json"
                    },
                    {
                        "url": "http://localhost:8000/cell-sets",
                        "fileType": "cell-sets.json"
                    }
                ]
            },
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_set_layout_single_view():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')
    my_view = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    vc.layout(my_view)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 0,
                'y': 0,
                'h': 12,
                'w': 12,
            }
        ],
        "initStrategy": "auto"
    }


def test_config_set_layout_multi_view():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')
    v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc.layout(hconcat(v1, vconcat(v2, v3)))

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 0,
                'y': 0,
                'h': 12,
                'w': 6,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 6,
                'y': 0,
                'h': 6,
                'w': 6,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 6,
                'y': 6,
                'h': 6,
                'w': 6,
            }
        ],
        "initStrategy": "auto"
    }


def test_config_set_layout_multi_view_custom_split():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')
    v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc.layout(hconcat(v1, vconcat(v2, v3, split=[1, 2]), split=[3, 1]))

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 0,
                'y': 0,
                'h': 12,
                'w': 9,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 9,
                'y': 0,
                'h': 4,
                'w': 3,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 9,
                'y': 4,
                'h': 8,
                'w': 3,
            }
        ],
        "initStrategy": "auto"
    }


def test_config_set_layout_multi_view_magic():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')
    v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc.layout(v1 | (v2 / v3))

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
        },
        "layout": [
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 0,
                'y': 0,
                'h': 12,
                'w': 6,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 6,
                'y': 0,
                'h': 6,
                'w': 6,
            },
            {
                'component': 'spatial',
                'coordinationScopes': {
                    'dataset': 'A',
                },
                'x': 6,
                'y': 6,
                'h': 6,
                'w': 6,
            }
        ],
        "initStrategy": "auto"
    }


def test_config_from_dict():
    vc = VitessceConfig.from_dict({
        "version": "1.0.15",
        "name": "Test name",
        "description": "Test description",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My First Dataset',
                'files': [
                    {
                        'url': 'http://cells.json',
                        'fileType': 'cells.json',
                        'requestInit': {
                            'headers': {
                                'Authorization': 'Bearer token'
                            }
                        }
                    }
                ]
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A'
            },
            'spatialZoom': {
                'ABC': 11
            },
        },
        "layout": [
            {
                "component": "spatial",
                "props": {
                    "cellRadius": 50
                },
                "coordinationScopes": {
                    "spatialZoom": 'ABC'
                },
                "x": 5,
                "y": 0,
                "w": 4,
                "h": 4
            },
        ],
        "initStrategy": "auto"
    })

    vc.add_dataset(name='My Second Dataset')

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.15",
        "name": "Test name",
        "description": "Test description",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My First Dataset',
                'files': [
                    {
                        'url': 'http://cells.json',
                        'fileType': 'cells.json',
                        'requestInit': {
                            'headers': {
                                'Authorization': 'Bearer token'
                            }
                        }
                    }
                ]
            },
            {
                'uid': 'B',
                'name': 'My Second Dataset',
                'files': []
            }
        ],
        'coordinationSpace': {
            'dataset': {
                'A': 'A',
                'B': 'B'
            },
            'spatialZoom': {
                'ABC': 11,
            },
        },
        "layout": [
            {
                "component": "spatial",
                "props": {
                    "cellRadius": 50
                },
                "coordinationScopes": {
                    "spatialZoom": 'ABC'
                },
                "x": 5,
                "y": 0,
                "w": 4,
                "h": 4
            },
        ],
        "initStrategy": "auto"
    }


def test_config_from_dict_raises_error_if_dataset_ambiguous():
    with pytest.raises(ValueError):
        VitessceConfig.from_dict({
            "version": "1.0.15",
            "name": "Test name",
            "description": "Test description",
            "datasets": [
                {
                    'uid': 'A',
                    'name': 'My First Dataset',
                    'files': [
                        {
                            'url': 'http://cells-1.json',
                            'fileType': 'cells.json'
                        }
                    ]
                },
                {
                    'uid': 'B',
                    'name': 'My Second Dataset',
                    'files': [
                        {
                            'url': 'http://cells-2.json',
                            'fileType': 'cells.json'
                        }
                    ]
                }
            ],
            'coordinationSpace': {
                'dataset': {
                    'A': 'A'
                },
                'spatialZoom': {
                    'ABC': 11
                },
            },
            "layout": [
                {
                    "component": "spatial",
                    "props": {
                        "cellRadius": 50
                    },
                    "coordinationScopes": {
                        "spatialZoom": 'ABC'
                    },
                    "x": 5,
                    "y": 0,
                    "w": 4,
                    "h": 4
                },
            ],
            "initStrategy": "auto"
        })


def test_config_to_python_with_data_objects():
    vc = VitessceConfig(schema_version="1.0.15")

    class MockWrapperA(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self._repr = make_repr(locals())
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i, base_dir=None):
            def get_molecules(base_url):
                return {
                    "url": f"{base_url}/molecules",
                    "fileType": "molecules.json",
                    "coordinationValues": {
                        "obsType": "molecule"
                    }
                }

            def get_cells(base_url):
                return {
                    "url": f"{base_url}/cells",
                    "fileType": "cells.json"
                }
            self.file_def_creators += [get_molecules, get_cells]

    class MockWrapperB(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self._repr = make_repr(locals())
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i, base_dir=None):
            def get_cell_sets(base_url):
                return {
                    "url": f"{base_url}/cell-sets",
                    "fileType": "cell-sets.json"
                }
            self.file_def_creators += [get_cell_sets]

    dataset_a = vc.add_dataset(name='My First Dataset').add_object(
        obj=MockWrapperA("Experiment A")
    ).add_file(
        url="http://example.com/my_cells.json",
        file_type=ft.CELLS_JSON
    )
    dataset_b = vc.add_dataset(name='My Second Dataset').add_object(
        obj=MockWrapperB("Experiment B")
    )
    vc.add_view(cm.SPATIAL, dataset=dataset_a, x=0, y=0,
                w=3, h=3).set_props(title="My spatial plot")
    vc.add_view(cm.SCATTERPLOT, dataset=dataset_b, x=3, y=0, w=3,
                h=3, mapping="PCA").set_props(title="My scatterplot")
    base_url = "http://localhost:8000"

    classes_to_import, code_block = vc.to_python()
    assert classes_to_import == ['VitessceChainableConfig', 'VitessceConfigDatasetFile']

    # Evaluate the code string directly
    reconstructed_vc = eval(code_block)
    assert vc.to_dict(base_url=base_url) == reconstructed_vc.to_dict(base_url=base_url)

    # Convert code string to an AST and back before evaluation
    if hasattr(ast, 'unparse'):  # pragma: no cover
        # Unparse added in Python 3.9
        ast_reconstructed_vc = eval(ast.unparse(ast.parse(code_block)))
        assert vc.to_dict(base_url=base_url) == ast_reconstructed_vc.to_dict(base_url=base_url)
    else:  # pragma: no cover
        ast.parse(code_block)


def test_config_add_coordination_by_dict():
    vc = VitessceConfig(schema_version="1.0.16")

    vc.add_coordination_by_dict({
        'spatialImageLayer': CL([
            {
                'image': 'S-1905-017737_bf',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialImageChannel': CL([
                    {
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'spatialTargetC': 1,
                        'spatialChannelColor': [0, 255, 0],
                    },
                ]),
            },
        ]),
        'spatialSegmentationLayer': CL([
            {
                'image': 'S-1905-017737',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialSegmentationChannel': CL([
                    {
                        'obsType': 'Cortical Interstitia',
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Non-Globally Sclerotic Glomeruli',
                        'spatialTargetC': 1,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Globally Sclerotic Glomeruli',
                        'spatialTargetC': 2,
                        'spatialChannelColor': [255, 0, 0],
                    },
                ]),
            },
        ]),
    })

    vc_dict = vc.to_dict()

    assert vc_dict["coordinationSpace"] == {
        "spatialImageLayer": {"A": '__dummy__'},
        "image": {"A": 'S-1905-017737_bf', "B": 'S-1905-017737'},
        "spatialLayerVisible": {"A": True, "B": True},
        "spatialLayerOpacity": {"A": 1, "B": 1},
        "spatialImageChannel": {"A": '__dummy__', "B": '__dummy__'},
        "spatialTargetC": {
            "A": 0, "B": 1, "C": 0, "D": 1, "E": 2,
        },
        "spatialChannelColor": {
            "A": [255, 0, 0], "B": [0, 255, 0], "C": [255, 0, 0], "D": [255, 0, 0], "E": [255, 0, 0],
        },
        "spatialSegmentationLayer": {"A": '__dummy__'},
        "spatialSegmentationChannel": {"A": '__dummy__', "B": '__dummy__', "C": '__dummy__'},
        "obsType": {
            "A": 'Cortical Interstitia',
            "B": 'Non-Globally Sclerotic Glomeruli',
            "C": 'Globally Sclerotic Glomeruli',
        },
    }


def test_config_add_and_use_coordination_by_dict():
    vc = VitessceConfig(schema_version="1.0.16", name="My config")
    dataset = vc.add_dataset(name="My dataset")

    (color_scope, ) = vc.add_coordination('spatialChannelColor')
    color_scope.set_value([255, 0, 0])

    scopes = vc.add_coordination_by_dict({
        "spatialImageLayer": CL([
            {
                "image": 'S-1905-017737_bf',
                "spatialLayerVisible": True,
                "spatialLayerOpacity": 1,
                "spatialImageChannel": CL([
                    {
                        "spatialTargetC": 0,
                        "spatialChannelColor": [0, 255, 0],
                    },
                    {
                        "spatialTargetC": 1,
                        "spatialChannelColor": [0, 0, 255],
                    },
                ]),
            },
        ]),
        "spatialSegmentationLayer": CL([
            {
                "image": 'S-1905-017737',
                "spatialLayerVisible": True,
                "spatialLayerOpacity": 1,
                "spatialSegmentationChannel": CL([
                    {
                        "obsType": 'Cortical Interstitia',
                        "spatialTargetC": 0,
                        "spatialChannelColor": color_scope,
                    },
                    {
                        "obsType": 'Non-Globally Sclerotic Glomeruli',
                        "spatialTargetC": 1,
                        "spatialChannelColor": color_scope,
                    },
                    {
                        "obsType": 'Globally Sclerotic Glomeruli',
                        "spatialTargetC": 2,
                        "spatialChannelColor": color_scope,
                    },
                ]),
            },
        ]),
    })

    spatial_view = vc.add_view('spatial', dataset=dataset)
    spatial_view.use_coordination_by_dict(scopes)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.16",
        "name": "My config",
        "description": "",
        "datasets": [
            {
                "uid": "A",
                "name": "My dataset",
                "files": []
            }
        ],
        "coordinationSpace": {
            "dataset": {
                "A": "A"
            },
            "spatialImageLayer": {
                "A": "__dummy__"
            },
            "image": {
                "A": "S-1905-017737_bf",
                "B": "S-1905-017737"
            },
            "spatialLayerVisible": {
                "A": True,
                "B": True
            },
            "spatialLayerOpacity": {
                "A": 1,
                "B": 1
            },
            "spatialImageChannel": {
                "A": "__dummy__",
                "B": "__dummy__"
            },
            "spatialTargetC": {
                "A": 0,
                "B": 1,
                "C": 0,
                "D": 1,
                "E": 2
            },
            "spatialChannelColor": {
                "A": [255, 0, 0],
                "B": [0, 255, 0],
                "C": [0, 0, 255]
            },
            "spatialSegmentationLayer": {
                "A": "__dummy__"
            },
            "spatialSegmentationChannel": {
                "A": "__dummy__",
                "B": "__dummy__",
                "C": "__dummy__"
            },
            "obsType": {
                "A": "Cortical Interstitia",
                "B": "Non-Globally Sclerotic Glomeruli",
                "C": "Globally Sclerotic Glomeruli"
            }
        },
        "layout": [
            {
                "component": "spatial",
                "coordinationScopes": {
                    "spatialImageLayer": ["A"],
                    "spatialSegmentationLayer": ["A"]
                },
                "coordinationScopesBy": {
                    "spatialImageLayer": {
                        "image": {
                            "A": "A"
                        },
                        "spatialLayerVisible": {
                            "A": "A"
                        },
                        "spatialLayerOpacity": {
                            "A": "A"
                        },
                        "spatialImageChannel": {
                            "A": ["A", "B"]
                        }
                    },
                    "spatialImageChannel": {
                        "spatialTargetC": {
                            "A": "A",
                            "B": "B"
                        },
                        "spatialChannelColor": {
                            "A": "B",
                            "B": "C"
                        }
                    },
                    "spatialSegmentationLayer": {
                        "image": {
                            "A": "B"
                        },
                        "spatialLayerVisible": {
                            "A": "B"
                        },
                        "spatialLayerOpacity": {
                            "A": "B"
                        },
                        "spatialSegmentationChannel": {
                            "A": ["A", "B", "C"]
                        }
                    },
                    "spatialSegmentationChannel": {
                        "obsType": {
                            "A": "A",
                            "B": "B",
                            "C": "C"
                        },
                        "spatialTargetC": {
                            "A": "C",
                            "B": "D",
                            "C": "E"
                        },
                        "spatialChannelColor": {
                            "A": "A",
                            "B": "A",
                            "C": "A"
                        }
                    }
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            }
        ],
        "initStrategy": "auto"
    }


def test_config_use_meta_complex_coordination():
    vc = VitessceConfig(schema_version="1.0.16", name="My config")
    dataset = vc.add_dataset(name="My dataset")

    scopes = vc.add_coordination_by_dict({
        'spatialImageLayer': CL([
            {
                'image': 'S-1905-017737_bf',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialImageChannel': CL([
                    {
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'spatialTargetC': 1,
                        'spatialChannelColor': [0, 255, 0],
                    },
                ]),
            },
        ]),
        'spatialSegmentationLayer': CL([
            {
                'image': 'S-1905-017737',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialSegmentationChannel': CL([
                    {
                        'obsType': 'Cortical Interstitia',
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Non-Globally Sclerotic Glomeruli',
                        'spatialTargetC': 1,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Globally Sclerotic Glomeruli',
                        'spatialTargetC': 2,
                        'spatialChannelColor': [255, 0, 0],
                    },
                ]),
            },
        ]),
    })

    meta_coordination_scope = vc.add_meta_coordination()
    meta_coordination_scope.use_coordination_by_dict(scopes)

    spatial_view = vc.add_view('spatial', dataset=dataset)
    lc_view = vc.add_view('layerController', dataset=dataset)

    spatial_view.use_meta_coordination(meta_coordination_scope)
    lc_view.use_meta_coordination(meta_coordination_scope)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.16",
        "name": "My config",
        "description": "",
        "datasets": [
            {
                "uid": "A",
                "name": "My dataset",
                "files": []
            }
        ],
        "coordinationSpace": {
            "dataset": {
                "A": "A"
            },
            "spatialImageLayer": {
                "A": "__dummy__"
            },
            "image": {
                "A": "S-1905-017737_bf",
                "B": "S-1905-017737"
            },
            "spatialLayerVisible": {
                "A": True,
                "B": True
            },
            "spatialLayerOpacity": {
                "A": 1,
                "B": 1
            },
            "spatialImageChannel": {
                "A": "__dummy__",
                "B": "__dummy__"
            },
            "spatialTargetC": {
                "A": 0,
                "B": 1,
                "C": 0,
                "D": 1,
                "E": 2
            },
            "spatialChannelColor": {
                "A": [255, 0, 0],
                "B": [0, 255, 0],
                "C": [255, 0, 0],
                "D": [255, 0, 0],
                "E": [255, 0, 0]
            },
            "spatialSegmentationLayer": {
                "A": "__dummy__"
            },
            "spatialSegmentationChannel": {
                "A": "__dummy__",
                "B": "__dummy__",
                "C": "__dummy__"
            },
            "obsType": {
                "A": "Cortical Interstitia",
                "B": "Non-Globally Sclerotic Glomeruli",
                "C": "Globally Sclerotic Glomeruli"
            },
            "metaCoordinationScopes": {
                "A": {
                    "spatialImageLayer": ["A"],
                    "spatialSegmentationLayer": ["A"]
                }
            },
            "metaCoordinationScopesBy": {
                "A": {
                    "spatialImageLayer": {
                        "image": {
                            "A": "A"
                        },
                        "spatialLayerVisible": {
                            "A": "A"
                        },
                        "spatialLayerOpacity": {
                            "A": "A"
                        },
                        "spatialImageChannel": {
                            "A": ["A", "B"]
                        }
                    },
                    "spatialImageChannel": {
                        "spatialTargetC": {
                            "A": "A",
                            "B": "B"
                        },
                        "spatialChannelColor": {
                            "A": "A",
                            "B": "B"
                        }
                    },
                    "spatialSegmentationLayer": {
                        "image": {
                            "A": "B"
                        },
                        "spatialLayerVisible": {
                            "A": "B"
                        },
                        "spatialLayerOpacity": {
                            "A": "B"
                        },
                        "spatialSegmentationChannel": {
                            "A": ["A", "B", "C"]
                        }
                    },
                    "spatialSegmentationChannel": {
                        "obsType": {
                            "A": "A",
                            "B": "B",
                            "C": "C"
                        },
                        "spatialTargetC": {
                            "A": "C",
                            "B": "D",
                            "C": "E"
                        },
                        "spatialChannelColor": {
                            "A": "C",
                            "B": "D",
                            "C": "E"
                        }
                    }
                }
            }
        },
        "layout": [
            {
                "component": "spatial",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["A"],
                    "metaCoordinationScopesBy": ["A"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            },
            {
                "component": "layerController",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["A"],
                    "metaCoordinationScopesBy": ["A"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            }
        ],
        "initStrategy": "auto"
    }


def test_config_link_views_by_dict():
    vc = VitessceConfig(schema_version="1.0.16", name="My config")
    dataset = vc.add_dataset(name="My dataset")

    spatial_view = vc.add_view('spatial', dataset=dataset)
    lc_view = vc.add_view('layerController', dataset=dataset)

    vc.link_views_by_dict([spatial_view, lc_view], {
        'spatialImageLayer': CL([
            {
                'image': 'S-1905-017737_bf',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialImageChannel': CL([
                    {
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'spatialTargetC': 1,
                        'spatialChannelColor': [0, 255, 0],
                    },
                ]),
            },
        ]),
        'spatialSegmentationLayer': CL([
            {
                'image': 'S-1905-017737',
                'spatialLayerVisible': True,
                'spatialLayerOpacity': 1,
                'spatialSegmentationChannel': CL([
                    {
                        'obsType': 'Cortical Interstitia',
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Non-Globally Sclerotic Glomeruli',
                        'spatialTargetC': 1,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'obsType': 'Globally Sclerotic Glomeruli',
                        'spatialTargetC': 2,
                        'spatialChannelColor': [255, 0, 0],
                    },
                ]),
            },
        ]),
    })

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.16",
        "name": "My config",
        "description": "",
        "datasets": [
            {
                "uid": "A",
                "name": "My dataset",
                "files": []
            }
        ],
        "coordinationSpace": {
            "dataset": {
                "A": "A"
            },
            "spatialImageLayer": {
                "A": "__dummy__"
            },
            "image": {
                "A": "S-1905-017737_bf",
                "B": "S-1905-017737"
            },
            "spatialLayerVisible": {
                "A": True,
                "B": True
            },
            "spatialLayerOpacity": {
                "A": 1,
                "B": 1
            },
            "spatialImageChannel": {
                "A": "__dummy__",
                "B": "__dummy__"
            },
            "spatialTargetC": {
                "A": 0,
                "B": 1,
                "C": 0,
                "D": 1,
                "E": 2
            },
            "spatialChannelColor": {
                "A": [255, 0, 0],
                "B": [0, 255, 0],
                "C": [255, 0, 0],
                "D": [255, 0, 0],
                "E": [255, 0, 0]
            },
            "spatialSegmentationLayer": {
                "A": "__dummy__"
            },
            "spatialSegmentationChannel": {
                "A": "__dummy__",
                "B": "__dummy__",
                "C": "__dummy__"
            },
            "obsType": {
                "A": "Cortical Interstitia",
                "B": "Non-Globally Sclerotic Glomeruli",
                "C": "Globally Sclerotic Glomeruli"
            },
            "metaCoordinationScopes": {
                "A": {
                    "spatialImageLayer": ["A"],
                    "spatialSegmentationLayer": ["A"]
                }
            },
            "metaCoordinationScopesBy": {
                "A": {
                    "spatialImageLayer": {
                        "image": {
                            "A": "A"
                        },
                        "spatialLayerVisible": {
                            "A": "A"
                        },
                        "spatialLayerOpacity": {
                            "A": "A"
                        },
                        "spatialImageChannel": {
                            "A": ["A", "B"]
                        }
                    },
                    "spatialImageChannel": {
                        "spatialTargetC": {
                            "A": "A",
                            "B": "B"
                        },
                        "spatialChannelColor": {
                            "A": "A",
                            "B": "B"
                        }
                    },
                    "spatialSegmentationLayer": {
                        "image": {
                            "A": "B"
                        },
                        "spatialLayerVisible": {
                            "A": "B"
                        },
                        "spatialLayerOpacity": {
                            "A": "B"
                        },
                        "spatialSegmentationChannel": {
                            "A": ["A", "B", "C"]
                        }
                    },
                    "spatialSegmentationChannel": {
                        "obsType": {
                            "A": "A",
                            "B": "B",
                            "C": "C"
                        },
                        "spatialTargetC": {
                            "A": "C",
                            "B": "D",
                            "C": "E"
                        },
                        "spatialChannelColor": {
                            "A": "C",
                            "B": "D",
                            "C": "E"
                        }
                    }
                }
            }
        },
        "layout": [
            {
                "component": "spatial",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["A"],
                    "metaCoordinationScopesBy": ["A"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            },
            {
                "component": "layerController",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["A"],
                    "metaCoordinationScopesBy": ["A"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            }
        ],
        "initStrategy": "auto"
    }


def test_config_link_views_by_dict_with_scope_prefix():
    vc = VitessceConfig(schema_version="1.0.16", name="My config")
    dataset = vc.add_dataset(name="My dataset")

    spatial_view = vc.add_view('spatial', dataset=dataset)
    lc_view = vc.add_view('layerController', dataset=dataset)

    vc.link_views_by_dict([spatial_view, lc_view], {
        'spatialImageLayer': CL([
            {
                'spatialLayerOpacity': 1,
                'spatialImageChannel': CL([
                    {
                        'spatialTargetC': 0,
                        'spatialChannelColor': [255, 0, 0],
                    },
                    {
                        'spatialTargetC': 1,
                        'spatialChannelColor': [0, 255, 0],
                    },
                ]),
            },
        ])
    }, scope_prefix="SOME_PREFIX_")

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.16",
        "name": "My config",
        "description": "",
        "datasets": [
            {
                "uid": "A",
                "name": "My dataset",
                "files": []
            }
        ],
        "coordinationSpace": {
            "dataset": {
                "A": "A"
            },
            "spatialImageLayer": {
                "SOME_PREFIX_0": "__dummy__"
            },
            "spatialLayerOpacity": {
                "SOME_PREFIX_0": 1,
            },
            "spatialImageChannel": {
                "SOME_PREFIX_0": "__dummy__",
                "SOME_PREFIX_1": "__dummy__"
            },
            "spatialTargetC": {
                "SOME_PREFIX_0": 0,
                "SOME_PREFIX_1": 1
            },
            "spatialChannelColor": {
                "SOME_PREFIX_0": [255, 0, 0],
                "SOME_PREFIX_1": [0, 255, 0],
            },
            "metaCoordinationScopes": {
                "SOME_PREFIX_0": {
                    "spatialImageLayer": ["SOME_PREFIX_0"]
                }
            },
            "metaCoordinationScopesBy": {
                "SOME_PREFIX_0": {
                    "spatialImageLayer": {
                        "spatialLayerOpacity": {
                            "SOME_PREFIX_0": "SOME_PREFIX_0"
                        },
                        "spatialImageChannel": {
                            "SOME_PREFIX_0": ["SOME_PREFIX_0", "SOME_PREFIX_1"]
                        }
                    },
                    "spatialImageChannel": {
                        "spatialTargetC": {
                            "SOME_PREFIX_0": "SOME_PREFIX_0",
                            "SOME_PREFIX_1": "SOME_PREFIX_1"
                        },
                        "spatialChannelColor": {
                            "SOME_PREFIX_0": "SOME_PREFIX_0",
                            "SOME_PREFIX_1": "SOME_PREFIX_1"
                        }
                    }
                }
            }
        },
        "layout": [
            {
                "component": "spatial",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["SOME_PREFIX_0"],
                    "metaCoordinationScopesBy": ["SOME_PREFIX_0"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            },
            {
                "component": "layerController",
                "coordinationScopes": {
                    "dataset": "A",
                    "metaCoordinationScopes": ["SOME_PREFIX_0"],
                    "metaCoordinationScopesBy": ["SOME_PREFIX_0"]
                },
                "x": 0, "y": 0, "w": 1, "h": 1
            }
        ],
        "initStrategy": "auto"
    }
