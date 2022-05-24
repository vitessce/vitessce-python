import pytest
import ast

from vitessce import (  # noqa: F401
    VitessceConfig,
    CoordinationType as ct,
    Component as cm,
    DataType as dt,
    FileType as ft,
    hconcat,
    vconcat,
    AbstractWrapper,
    make_repr,

    # Neither of these is in the source code, but they do appear in code which is eval'd.
    VitessceChainableConfig,
    VitessceConfigDatasetFile,
)


def test_config_creation():
    vc = VitessceConfig()
    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
        "name": "",
        "description": "",
        "datasets": [],
        "coordinationSpace": {},
        "layout": [],
        "initStrategy": "auto"
    }


def test_config_add_dataset():
    vc = VitessceConfig()
    vc.add_dataset(name='My Dataset')

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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


def test_config_add_dataset_add_files():
    vc = VitessceConfig()
    vc.add_dataset(name='My Chained Dataset').add_file(
        url="http://example.com/cells.json",
        data_type=dt.CELLS,
        file_type=ft.CELLS_JSON,
    ).add_file(
        url="http://example.com/cell_sets.json",
        data_type=dt.CELL_SETS,
        file_type=ft.CELL_SETS_JSON,
    )

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Chained Dataset',
                'files': [
                    {
                        'url': 'http://example.com/cells.json',
                        'type': 'cells',
                        'fileType': 'cells.json'
                    },
                    {
                        'url': 'http://example.com/cell_sets.json',
                        'type': 'cell-sets',
                        'fileType': 'cell-sets.json'
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
    vc = VitessceConfig()
    my_dataset = vc.add_dataset(name='My Dataset')

    vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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
    vc = VitessceConfig()
    my_dataset = vc.add_dataset(name='My Dataset')

    vc.add_view(
        cm.SCATTERPLOT, dataset=my_dataset, mapping="X_umap")

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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
    vc = VitessceConfig()
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
        "version": "1.0.7",
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
    vc = VitessceConfig()

    class MockWrapperA(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i):
            def get_molecules(base_url):
                return {
                    "url": f"{base_url}/molecules",
                    "type": "molecules",
                    "fileType": "molecules.json"
                }

            def get_cells(base_url):
                return {
                    "url": f"{base_url}/cells",
                    "type": "cells",
                    "fileType": "cells.json"
                }
            self.file_def_creators += [get_molecules, get_cells]

    class MockWrapperB(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i):
            def get_cell_sets(base_url):
                return {
                    "url": f"{base_url}/cell-sets",
                    "type": "cell-sets",
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
        "version": "1.0.7",
        "name": "",
        "description": "",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My Object Dataset',
                'files': [
                    {
                        "url": "http://localhost:8000/molecules",
                        "type": "molecules",
                        "fileType": "molecules.json"
                    },
                    {
                        "url": "http://localhost:8000/cells",
                        "type": "cells",
                        "fileType": "cells.json"
                    },
                    {
                        "url": "http://localhost:8000/cell-sets",
                        "type": "cell-sets",
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
    vc = VitessceConfig()
    my_dataset = vc.add_dataset(name='My Dataset')
    my_view = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    vc.layout(my_view)

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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
    vc = VitessceConfig()
    my_dataset = vc.add_dataset(name='My Dataset')
    v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc.layout(hconcat(v1, vconcat(v2, v3)))

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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


def test_config_set_layout_multi_view_magic():
    vc = VitessceConfig()
    my_dataset = vc.add_dataset(name='My Dataset')
    v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
    v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)

    vc.layout(v1 | (v2 / v3))

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
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
        "version": "1.0.7",
        "name": "Test name",
        "description": "Test description",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My First Dataset',
                'files': [
                    {
                        'url': 'http://cells.json',
                        'type': 'cells',
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

    vc.add_dataset(name='My Second Dataset')

    vc_dict = vc.to_dict()

    assert vc_dict == {
        "version": "1.0.7",
        "name": "Test name",
        "description": "Test description",
        "datasets": [
            {
                'uid': 'A',
                'name': 'My First Dataset',
                'files': [
                    {
                        'url': 'http://cells.json',
                        'type': 'cells',
                        'fileType': 'cells.json'
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
            "version": "1.0.7",
            "name": "Test name",
            "description": "Test description",
            "datasets": [
                {
                    'uid': 'A',
                    'name': 'My First Dataset',
                    'files': [
                        {
                            'url': 'http://cells-1.json',
                            'type': 'cells',
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
                            'type': 'cells',
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
    vc = VitessceConfig()

    class MockWrapperA(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self._repr = make_repr(locals())
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i):
            def get_molecules(base_url):
                return {
                    "url": f"{base_url}/molecules",
                    "type": "molecules",
                    "fileType": "molecules.json"
                }

            def get_cells(base_url):
                return {
                    "url": f"{base_url}/cells",
                    "type": "cells",
                    "fileType": "cells.json"
                }
            self.file_def_creators += [get_molecules, get_cells]

    class MockWrapperB(AbstractWrapper):
        def __init__(self, name, **kwargs):
            super().__init__(**kwargs)
            self._repr = make_repr(locals())
            self.name = name

        def convert_and_save(self, dataset_uid, obj_i):
            def get_cell_sets(base_url):
                return {
                    "url": f"{base_url}/cell-sets",
                    "type": "cell-sets",
                    "fileType": "cell-sets.json"
                }
            self.file_def_creators += [get_cell_sets]

    dataset_a = vc.add_dataset(name='My First Dataset').add_object(
        obj=MockWrapperA("Experiment A")
    ).add_file(
        url="http://example.com/my_cells.json",
        file_type=ft.CELLS_JSON,
        data_type=dt.CELLS
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
