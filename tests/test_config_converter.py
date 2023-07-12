import pytest
from unittest.mock import patch, Mock
import os
from copy import deepcopy

from vitessce import (CellBrowserToVitessceConfigConverter, convert)
from vitessce.data_utils import (
    VAR_CHUNK_SIZE
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
      "SPINK1",
      "TRY6",
      "PNLIPRP1",
      "PLA2G1B",
      "CTRB1"
    ],
    "alpha": [
      "CRYBA2",
      "FAP",
      "APOH",
      "PALLD",
      "ST18"
    ],
    "beta": [
      "IAPP",
      "UCHL1",
      "HADH",
      "ERO1LB",
      "ABCC8"
    ],
    "delta": [
      "SST",
      "RBP4",
      "LEPR",
      "HADH",
      "RGS2"
    ],
    "ductal": [
      "ANXA4",
      "CLDN1",
      "SERPING1",
      "CFTR",
      "S100A10"
    ],
    "mesenchymal": [
      "SPARCL1",
      "SPARC",
      "MGP",
      "THBS1",
      "TIMP3"
    ],
    "nan": [
      "ERCC-00092",
      "ERCC-00022",
      "ERCC-00076",
      "ERCC-00112",
      "ERCC-00131"
    ],
    "unsure": [
      "G6PC2",
      "PCSK1",
      "CASR",
      "ELMO1",
      "GAD2"
    ]
  },
}

invalid_cellbrowser_config = deepcopy(valid_cellbrowser_config["fileVersions"])

@pytest.fixture
def mock_requests_get():
    with patch('requests.get') as mock_get:
        
        yield mock_get

def test_download_valid_config(mock_requests_get):

    # Set up the Mock to return a fake response when called
    mock_response = Mock()
    mock_response.json.return_value = valid_cellbrowser_config
    mock_requests_get.return_value = mock_response
    
    obj = CellBrowserToVitessceConfigConverter("test-project", "test-output-dir", False)
    is_valid = obj.download_config()

    # Now you can make assertions about how the mock was used and the result of your function
    mock_requests_get.assert_called_once_with('https://cells.ucsc.edu/test-project/dataset.json')
    assert is_valid == True
    assert obj.cellbrowser_config == valid_cellbrowser_config


def test_download_invalid_config(mock_requests_get):
  
    # Set up the Mock to return a fake response when called
    mock_response = Mock()
    mock_response.json.return_value = invalid_cellbrowser_config
    mock_requests_get.return_value = mock_response
    
    obj = CellBrowserToVitessceConfigConverter("test-project", "test-output-dir", False)
    is_valid = obj.download_config()

    # Now you can make assertions about how the mock was used and the result of your function
    mock_requests_get.assert_called_once_with('https://cells.ucsc.edu/test-project/dataset.json')
    assert is_valid == False
    assert obj.cellbrowser_config == invalid_cellbrowser_config


@pytest.fixture
def mock_end_to_end_tests():
    # Set up the Mock to return a fake response when called
    mock_response_json = Mock()
    mock_response_json.json.return_value = valid_cellbrowser_config
    mock_response_json.raise_for_status.return_value = None
    mock_response_json.content = b''

    with open('tests/data/smaller_expr_matrix.tsv.gz', 'rb') as f:
        mock_response_expr_matrix = Mock()
        mock_response_expr_matrix.content = f.read()
        mock_response_expr_matrix.raise_for_status.return_value = None

    with open('tests/data/test_meta.tsv', 'rb') as f:
        mock_response_meta = Mock()
        mock_response_meta.content = f.read()
        mock_response_meta.raise_for_status.return_value = None

    with open('tests/data/test.coords.tsv.gz', 'rb') as f:
        mock_response_coords = Mock()
        mock_response_coords.content = f.read()
        mock_response_coords.raise_for_status.return_value = None

    with patch('requests.get') as mock_get:
        mock_get.side_effect = [mock_response_json, mock_response_expr_matrix, mock_response_meta, mock_response_coords]
        yield mock_get


@pytest.fixture
def mock_filter_cells():
    with patch('scanpy.pp.filter_cells') as mock:
        yield mock

@pytest.fixture
def mock_makedirs():
    with patch('os.makedirs') as mock:
        yield mock

# Define a fixture for adata.write_zarr
@pytest.fixture
def mock_write_zarr():
    with patch('anndata.AnnData.write_zarr') as mock:
        yield mock

def test_end_to_end(mock_makedirs, mock_write_zarr, mock_filter_cells, mock_end_to_end_tests):
    convert("test-project", "iva", keep_only_marker_genes=False)

    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/dataset.json")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/exprMatrix.tsv.gz")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/meta.tsv")
    mock_end_to_end_tests.assert_any_call("https://cells.ucsc.edu/test-project/test.coords.tsv.gz")

    # Assert the number of times requests.get was called
    assert mock_end_to_end_tests.call_count == 4
    assert mock_filter_cells.call_count == 1
    mock_makedirs.assert_called_once_with(os.path.dirname("iva/test-project"), exist_ok=True)
    mock_write_zarr.assert_called_once_with("iva/test-project/out.adata.zarr", chunks=[8, VAR_CHUNK_SIZE])
