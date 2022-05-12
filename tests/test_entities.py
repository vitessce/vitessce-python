import pytest

from vitessce import (
    CellSets,
    Cells,
)


def test_cells():

    cell_ids = ['cell_1', 'cell_2', 'cell_3']

    cells = Cells(cell_ids=cell_ids)
    assert list(cells.json.keys()) == cell_ids

    cells.add_mapping('umap', [[1, 1], [2, 2], [3, 3]])

    cells.add_mapping('pca', [[1, 1], [2, 2], [3, 3]])

    cells.add_centroids([[1, 1], [2, 2], [3, 3]])

    cells.add_polygon_outline(
        [[[1, 1], [1, 1], [1, 1]], [[2, 2], [2, 2], [2, 2]], [[3, 3], [3, 3], [3, 3]]])

    assert cells.json == {
        'cell_1': {
            'mappings': {'umap': [1, 1], 'pca': [1, 1]},
            'xy': [1, 1],
            'poly': [[1, 1], [1, 1], [1, 1]]
        },
        'cell_2': {
            'mappings': {'umap': [2, 2], 'pca': [2, 2]},
            'xy': [2, 2],
            'poly': [[2, 2], [2, 2], [2, 2]]
        },
        'cell_3': {
            'mappings': {'umap': [3, 3], 'pca': [3, 3]},
            'xy': [3, 3],
            'poly': [[3, 3], [3, 3], [3, 3]]
        }
    }


def test_cells_bad_polygon_outline_type():

    cell_ids = ['cell_1', 'cell_2', 'cell_3']
    cells = Cells(cell_ids=cell_ids)
    with pytest.raises(Exception) as context:
        # The extra 3 should be problematic since polygons are two dimensional.
        cells.add_polygon_outline([
            [[1, 1, 3], [1, 1], [1, 1]],
            [[2, 2], [2, 2], [2, 2]],
            [[3, 3], [3, 3], [3, 3]]
        ])
    assert 'Polygon outline for cell_1 should be a list of two element lists i.e xy coordinates' in str(context)


def test_cells_bad_mappings_length():

    cell_ids = ['cell_1', 'cell_2', 'cell_3']
    cells = Cells(cell_ids=cell_ids)
    with pytest.raises(Exception) as context:
        # There are 3 cells in this object so only passing in two scatterplot cooridnates is problematic.
        cells.add_mapping('umap', [[1, 1], [2, 2]])
    assert 'Coordinates length does not match Cell IDs Length' in str(context)


def test_cell_sets():

    cell_sets = CellSets()
    cell_sets.add_level_zero_node('Clusters')

    cell_sets.add_node('Cluster 1', ['Clusters'])
    cell_sets.add_node('Cluster 2', ['Clusters'])
    cell_sets.add_node('Subcluster 1', ['Clusters', 'Cluster 1'], ['cell_1', 'cell_2'])
    cell_sets.add_node('Subcluster 2', ['Clusters', 'Cluster 1'], ['cell_3', 'cell_4'])

    assert cell_sets.json == {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": [{
            "name": 'Clusters',
            "children": [
                {
                    "name": 'Cluster 1',
                    "children": [
                        {
                            "name": 'Subcluster 1',
                            "set": ['cell_1', 'cell_2']
                        },
                        {
                            "name": 'Subcluster 2',
                            "set": ['cell_3', 'cell_4']
                        }
                    ]
                },
                {
                    "name": 'Cluster 2',
                }
            ]
        }]
    }


def test_cell_sets_node_not_found():

    cell_sets = CellSets()
    with pytest.raises(Exception) as context:
        cell_sets.add_node('Cluster 1', ['Clusters Not Found'])

    assert "No node with path ['Clusters Not Found'] found to add Cluster 1 to" in str(context)
