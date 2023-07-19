import pytest
from unittest.mock import patch, Mock
from copy import deepcopy

from vitessce import (
    CellBrowserToAnndataZarrConverter,
    convert_cell_browser_project_to_anndata,
)


valid_cellbrowser_config = {
    "fileVersions": {
        "inMatrix": {
            "fname": "/hive/data/inside/cells/datasets/adultPancreas/exprMatrix.csv.gz",
            "md5": "8da7a759a8",
            "size": 23363664,
            "mtime": "2018-10-16 23:29:40"
        },
        "outMatrix": {
            "fname": "/usr/local/apache/htdocs-cells/adultPancreas/exprMatrix.tsv.gz",
            "md5": "934bbdeacd",
            "size": 22710325,
            "mtime": "2022-05-18 22:34:06"
        },
        "inMeta": {
            "fname": "/hive/data/inside/cells/datasets/adultPancreas/meta.tsv",
            "md5": "7699cf188d",
            "size": 527639,
            "mtime": "2019-02-26 16:08:50"
        },
        "outMeta": {
            "fname": "/usr/local/apache/htdocs-cells/adultPancreas/meta.tsv",
            "md5": "cdfeda9e0a",
            "size": 522326,
            "mtime": "2022-05-24 18:01:35"
        },
    },
    "coords": [
        {
            "name": "coords_0",
            "shortLabel": "t-SNE",
            "md5": "3ff37334ef",
            "minX": 0,
            "maxX": 65535,
            "minY": 0,
            "maxY": 65535,
            "type": "Uint16",
            "textFname": "test.coords.tsv.gz",
            "labelMd5": "d41d8cd98f"
        }
    ],
    "topMarkers": {
        "acinar": [
            "A1CF",
        ],
        "alpha": [
            "A1BG-AS1",
        ],
        "beta": [
            "LEPR",
        ],
        "delta": [
            "SST",
            "RBP4",
        ],
        "ductal": [
            "ANXA4",
        ],
        "mesenchymal": [
            "SPARCL1",
        ],
        "nan": [
            "ERCC-00092",
        ],
        "unsure": [
            "G6PC2",
            "PCSK1",
        ]
    },
}

invalid_cellbrowser_config = deepcopy(valid_cellbrowser_config["fileVersions"])

coords_with_no_fname = {
    "coords": [
        {
            "name": "coords_0",
            "shortLabel": "t-SNE",
            "md5": "3ff37334ef",
            "minX": 0,
            "maxX": 65535,
            "minY": 0,
            "maxY": 65535,
            "type": "Uint16",
            "labelMd5": "d41d8cd98f"
        }
    ],
}

coords_with_no_fname_multi_word_shortlabel = {
    "coords": [
        {
            "name": "coords_0",
            "shortLabel": "Seurat umap",
            "md5": "3ff37334ef",
            "minX": 0,
            "maxX": 65535,
            "minY": 0,
            "maxY": 65535,
            "type": "Uint16",
            "labelMd5": "d41d8cd98f"
        }
    ],
}

cellbrowser_config_no_coords_filename = {
    **valid_cellbrowser_config,
    **coords_with_no_fname
}

cellbrowser_config_no_coords_filename_multi_word_shortlabel = {
    **valid_cellbrowser_config,
    **coords_with_no_fname_multi_word_shortlabel
}

project_name = "test-project"

"""
@pytest.fixture
def mock_makedirs():
    with patch('os.makedirs') as mock:
        yield mock


@pytest.fixture
def mock_write_zarr():
    with patch('anndata.AnnData.write_zarr') as mock:
        yield mock


@pytest.fixture
def mock_write_adata():
    with patch('anndata.AnnData.write') as mock:
        yield mock

@pytest.fixture
def mock_read_adata():
    with patch('anndata.read') as mock:
        # Create a dummy AnnData object
        dummy_data = np.random.rand(10, 10)
        mock_adata = anndata.AnnData(dummy_data)

        # Mock _read_adata to return the dummy AnnData object
        mock.return_value = mock_adata
        yield mock
"""


@pytest.fixture
def mock_filter_cells():
    with patch('scanpy.pp.filter_cells') as mock:
        yield mock


@pytest.fixture
def mock_end_to_end_tests(request):
    config = request.param
    mock_response_json = Mock()
    mock_response_json.json.return_value = config

    with open('tests/data/smaller_expr_matrix.tsv.gz', 'rb') as f:
        mock_response_expr_matrix = Mock()
        mock_response_expr_matrix.content = f.read()

    with open('tests/data/test_meta.tsv', 'rb') as f:
        mock_response_meta = Mock()
        mock_response_meta.content = f.read()

    with open('tests/data/test.coords.tsv.gz', 'rb') as f:
        mock_response_coords = Mock()
        mock_response_coords.content = f.read()

    with patch('requests.get') as mock_get:
        mock_get.side_effect = [mock_response_json, mock_response_expr_matrix, mock_response_meta, mock_response_coords]
        yield mock_get


def test_download_valid_config():

    with patch('requests.get') as mock_get:
        mock_get.return_value.json.return_value = valid_cellbrowser_config
        obj = CellBrowserToAnndataZarrConverter(project_name, keep_only_marker_genes=False)
        is_valid = obj.download_config()

        mock_get.assert_called_once_with('https://cells.ucsc.edu/test-project/dataset.json')
        assert is_valid
        assert obj.cellbrowser_config == valid_cellbrowser_config


@pytest.mark.parametrize('mock_end_to_end_tests', [valid_cellbrowser_config], indirect=True)
def test_filter_based_on_marker_genes(mock_end_to_end_tests, mock_filter_cells):

    inst = CellBrowserToAnndataZarrConverter(project_name, keep_only_marker_genes=True)
    config_is_valid = inst.download_config()

    assert config_is_valid
    inst.create_anndata_object()

    assert inst.adata.shape == (8, 1)

    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/meta.tsv")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/test.coords.tsv.gz")

    assert mock_end_to_end_tests.call_count == 4
    assert mock_filter_cells.call_count == 1


@pytest.mark.parametrize(
    'mock_end_to_end_tests, expected',
    [
        (valid_cellbrowser_config, "test.coords.tsv.gz"),
        (cellbrowser_config_no_coords_filename, "tMinusSNE.coords.tsv.gz"),
        (cellbrowser_config_no_coords_filename_multi_word_shortlabel, "Seurat_umap.coords.tsv.gz")
    ], indirect=["mock_end_to_end_tests"]
)
def test_end_to_end(mock_filter_cells, mock_end_to_end_tests, expected):
    convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/meta.tsv")
    mock_end_to_end_tests.assert_any_call(f"https://cells.ucsc.edu/test-project/{expected}")

    assert mock_end_to_end_tests.call_count == 4
    assert mock_filter_cells.call_count == 1


def test_end_to_end_invalid_config(mock_filter_cells):
    with patch('requests.get') as mock_get:
        mock_get.return_value.json.return_value = invalid_cellbrowser_config
        with pytest.raises(ValueError):
            convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

        mock_get.assert_called_once_with("https://cells.ucsc.edu/test-project/dataset.json")

    assert mock_get.call_count == 1
    assert mock_filter_cells.call_count == 0


def test_end_to_end_download_config_raises_exception(mock_filter_cells):
    mock_response = Mock()
    mock_response.raise_for_status.side_effect = Exception("Error downloading file")

    with patch('requests.get') as mock_get:
        mock_get.return_value = mock_response
        with pytest.raises(Exception):
            convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

        mock_get.assert_called_once_with("https://cells.ucsc.edu/test-project/dataset.json")

        assert mock_get.call_count == 1
    assert mock_filter_cells.call_count == 0


def test_end_to_end_load_expr_matrix_raises_exception(mock_filter_cells):
    mock_first_response = Mock()
    mock_first_response.json.return_value = valid_cellbrowser_config

    mock_second_response = Mock()
    mock_second_response.raise_for_status.side_effect = Exception("Error downloading file")

    with patch('requests.get') as mock_get:
        mock_get.side_effect = [mock_first_response, mock_second_response]
        with pytest.raises(Exception):
            convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
        assert mock_get.call_count == 2
    assert mock_filter_cells.call_count == 0


def test_end_to_end_load_cell_metadata_raises_exception(mock_filter_cells):
    mock_get_config = Mock()
    mock_get_config.json.return_value = valid_cellbrowser_config

    with open('tests/data/smaller_expr_matrix.tsv.gz', 'rb') as f:
        mock_response_expr_matrix = Mock()
        mock_response_expr_matrix.content = f.read()
        mock_response_expr_matrix.raise_for_status.return_value = None

    mock_response_meta = Mock()
    mock_response_meta.raise_for_status.side_effect = Exception("Error downloading file")
    assert mock_filter_cells.call_count == 0

    with patch('requests.get') as mock_get:
        mock_get.side_effect = [mock_get_config, mock_response_expr_matrix, mock_response_meta]
        with pytest.raises(Exception):
            convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/meta.tsv")
        assert mock_get.call_count == 3
    assert mock_filter_cells.call_count == 0


def test_end_to_end_add_coords_raises_exception(mock_filter_cells):
    mock_get_config = Mock()
    mock_get_config.json.return_value = valid_cellbrowser_config

    with open('tests/data/smaller_expr_matrix.tsv.gz', 'rb') as f:
        mock_response_expr_matrix = Mock()
        mock_response_expr_matrix.content = f.read()
        mock_response_expr_matrix.raise_for_status.return_value = None

    with open('tests/data/test_meta.tsv', 'rb') as f:
        mock_response_meta = Mock()
        mock_response_meta.content = f.read()
        mock_response_meta.raise_for_status.return_value = None

    mock_coords = Mock()
    mock_coords.raise_for_status.side_effect = Exception("Error downloading file")

    with patch('requests.get') as mock_get:
        mock_get.side_effect = [mock_get_config, mock_response_expr_matrix, mock_response_meta, mock_coords]
        with pytest.raises(Exception):
            convert_cell_browser_project_to_anndata(project_name, keep_only_marker_genes=False)

        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/meta.tsv")
        mock_get.assert_any_call("https://cells.ucsc.edu/test-project/test.coords.tsv.gz")
        assert mock_get.call_count == 4
    assert mock_filter_cells.call_count == 0
