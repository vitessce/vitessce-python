import json
from pathlib import Path

import pytest

from vitessce import (
    VitessceConfig,
    ViewType as vt,
    FileType as ft,
    CoordinationType as ct,
    validate_config,
)


@pytest.fixture
def config():
    vc = VitessceConfig(schema_version="1.0.18", name="My config")
    dataset = vc.add_dataset(name="My dataset")
    dataset.add_file(url="https://example.com/sets.json", file_type=ft.OBS_SETS_JSON)
    dataset.add_file(
        url="https://example.com/a.zarr",
        file_type=ft.OBS_EMBEDDING_ANNDATA_ZARR,
        options={"path": "obsm/X_umap"},
        coordination_values={"embeddingType": "UMAP"},
    )
    scatterplot = vc.add_view(vt.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    heatmap = vc.add_view(vt.HEATMAP, dataset=dataset)
    vc.link_views([scatterplot, heatmap], [ct.EMBEDDING_ZOOM], [2])
    vc.layout(scatterplot | heatmap)
    return json.loads(json.dumps(vc.to_dict(base_url="http://localhost:8000")))


def test_valid(config):
    assert validate_config(config) == []


def test_unknown_component(config):
    config["layout"][0]["component"] = "scatterplott"
    [error] = validate_config(config)
    assert error["path"] == "$.layout[0].component"
    assert error["message"].startswith("'scatterplott' is not one of ['description',")


def test_file_options_follow_file_type(config):
    config["datasets"][0]["files"][1]["options"] = {"pth": "obsm/X_umap"}
    assert validate_config(config) == [
        {"path": "$.datasets[0].files[1].options", "message": "'path' is a required property"},
        {"path": "$.datasets[0].files[1].options", "message": "Additional properties are not allowed ('pth' was unexpected)"},
    ]


def test_unknown_file_type(config):
    config["datasets"][0]["files"][1]["fileType"] = "obsEmbedding.anndata.zar"
    assert validate_config(config) == [
        {"path": "$.datasets[0].files[1].fileType", "message": "'obsEmbedding.anndata.zar' is not a known fileType"},
    ]


def test_coordination_value_skips_undefined_branch(config):
    config["coordinationSpace"]["embeddingZoom"]["A"] = "two"
    assert validate_config(config) == [
        {"path": "$.coordinationSpace.embeddingZoom.A", "message": "'two' is not of type 'number', 'null'"},
    ]


def test_nullable_coordination_value_follows_type():
    # Upgrading this fixture leaves 1.0.0-style spatial layers that the base plugins reject.
    fixtures = json.loads((Path(__file__).parent / "upgrade_fixtures.json").read_text())
    errors = validate_config(fixtures["legacyViewConfig1_0_0"])
    assert {"path": "$.coordinationSpace.spatialImageLayer.A[0]", "message": "'channels' is a required property"} in errors
    assert validate_config(fixtures["legacyViewConfig1_0_0"], base_plugins=False) == []


def test_ambiguous_any_of_keeps_original_error(config):
    config["layout"][0]["coordinationScopes"]["embeddingZoom"] = 5
    assert validate_config(config) == [
        {"path": "$.layout[0].coordinationScopes.embeddingZoom", "message": "5 is not valid under any of the given schemas"},
    ]


def test_older_version_is_upgraded(config):
    config["version"] = "1.0.15"
    assert validate_config(config) == []
    config["layout"][0]["component"] = "scatterplott"
    assert validate_config(config)[0]["path"] == "$.layout[0].component"


def test_older_version_structural_errors(config):
    config["version"] = "1.0.15"
    del config["layout"][0]["x"]
    assert validate_config(config) == [{"path": "$.layout[0]", "message": "'x' is a required property"}]


def test_upgrade_error(config):
    config["version"] = "1.0.17"
    config["datasets"][0]["files"].append({"fileType": "spatialdata.zarr", "options": {"labels": {}, "shapes": {}}})
    [error] = validate_config(config)
    assert error["path"] == "$"
    assert error["message"].startswith("Cannot upgrade both labels and shapes")


def test_plugin_types(config):
    config["layout"][0]["component"] = "myPluginView"
    assert validate_config(config) != []
    assert validate_config(config, base_plugins=False) == []


@pytest.mark.parametrize("bad_config, expected", [
    ({"version": "9.9.9"}, [{"path": "$.version", "message": "Config version was not recognized: '9.9.9'"}]),
    ({}, [{"path": "$.version", "message": "Config version was not recognized: None"}]),
    ([], [{"path": "$", "message": "config must be an object"}]),
])
def test_bad_input(bad_config, expected):
    assert validate_config(bad_config) == expected
