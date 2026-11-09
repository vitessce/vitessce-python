import json
from pathlib import Path

import pytest

from vitessce import upgrade_config
from vitessce.upgraders import (
    _upgrade_from_0_1_0,
    _upgrade_from_1_0_0,
    _upgrade_from_1_0_15,
    _upgrade_from_1_0_16,
)

# Exported from view-config-utils.test.fixtures.js in @vitessce/schemas@4.1.3.
FIXTURES = json.loads((Path(__file__).parent / "upgrade_fixtures.json").read_text())


def fixture(name):
    return json.loads(json.dumps(FIXTURES[name]))


def view(component, **kwargs):
    return {"component": component, "x": 0, "y": 0, **kwargs}


def files_config(version, *file_defs):
    return {"version": version, "name": "", "datasets": [{"uid": "A", "files": list(file_defs)}], "layout": []}


# Ports of the tests in view-config-utils.test.js.
# The JS upgradeFrom* functions also set the version; here upgrade_config does that.
@pytest.mark.parametrize("upgrade, before, after", [
    (lambda c: _upgrade_from_0_1_0(c, "A"), "legacyViewConfig0_1_0", "upgradedLegacyViewConfig0_1_0"),
    (_upgrade_from_1_0_0, "legacyViewConfig1_0_0", "upgradedLegacyViewConfig1_0_0"),
    (_upgrade_from_1_0_15, "implicitPerDatasetCoordinations", "explicitPerDatasetCoordinations"),
    (_upgrade_from_1_0_16, "annDataOptionsArrays", "annDataOptionsObjects"),
])
def test_upgrade_fixtures(upgrade, before, after):
    expected = fixture(after)
    assert {**upgrade(fixture(before)), "version": expected["version"]} == expected


def test_upgrades_more_than_once():
    config = fixture("legacyViewConfig1_0_0")
    upgraded = upgrade_config(config)
    assert upgraded["version"] == "1.0.18"
    assert config == FIXTURES["legacyViewConfig1_0_0"]


def test_upgrade_from_0_1_0_layer_controller_and_random_uid():
    config = {
        "version": "0.1.0",
        "layers": [],
        "staticLayout": [
            view("spatial", props={"view": {"zoom": 1, "target": [2, 3]}}),
            view("layerController"),
        ],
    }
    upgraded = upgrade_config(config)
    assert upgraded["layout"][1]["coordinationScopes"] == upgraded["layout"][0]["coordinationScopes"]
    assert upgraded["layout"][1]["coordinationScopes"] is not upgraded["layout"][0]["coordinationScopes"]
    assert len(upgraded["datasets"][0]["uid"]) == 36


def test_upgrade_from_1_0_1_and_1_0_2():
    config = {
        "version": "1.0.1",
        "datasets": [],
        "coordinationSpace": {"spatialRasterLayers": {"A": [{"type": "bitmask"}, {"type": "other"}], "B": None}},
        "layout": [view("layerController", props={"a": 1}), view("spatial")],
    }
    upgraded = upgrade_config(config)
    assert upgraded["coordinationSpace"]["spatialImageLayer"]["A"] == [{"type": "bitmask"}, {"type": "raster"}]
    assert upgraded["layout"][0]["props"] == {"a": 1, "globalDisable3d": True, "disableChannelsIfRgbDetected": True}
    assert "props" not in upgraded["layout"][1]


def test_upgrade_from_1_0_10():
    config = {
        "version": "1.0.10",
        "datasets": [],
        "coordinationSpace": {"cellFilter": {}, "geneSelection": None},
        "layout": [view("scatterplot", coordinationScopes={"cellHighlight": "A", "geneFilter": ""}), view("heatmap")],
    }
    upgraded = upgrade_config(config)
    # JS treats {} as truthy and null as falsy.
    assert upgraded["coordinationSpace"] == {"obsFilter": {}, "geneSelection": None}
    assert upgraded["layout"][0]["coordinationScopes"] == {"obsHighlight": "A", "geneFilter": ""}
    assert upgraded["layout"][1]["coordinationScopes"] == {}


def test_upgrade_from_1_0_11():
    config = files_config("1.0.11", {"fileType": "cells.json"}, {"fileType": "cell-sets.json"})
    config["coordinationSpace"] = {"embeddingType": {"A": "UMAP", "B": "PCA"}}
    files = upgrade_config(config)["datasets"][0]["files"]
    assert files[0]["options"] == {"embeddingTypes": ["UMAP", "PCA"]}
    assert "options" not in files[1]


def test_upgrade_from_1_0_12():
    config = {
        "version": "1.0.12",
        "coordinationSpace": {"dataset": {"A": "first", "B": "second"}, "obsLabelsType": {"A": "existing"}},
        "datasets": [
            {"uid": "first", "files": [{"fileType": "anndata-cells.zarr", "options": {"factors": ["obs/a", "obs/b"]}}]},
            {"uid": "second", "files": [{"fileType": "anndata-cells.zarr", "options": {}}]},
        ],
        "layout": [view("scatterplot", coordinationScopes={"dataset": "A"}), view("heatmap"), view("status", coordinationScopes={"dataset": "B"})],
    }
    upgraded = upgrade_config(config)
    assert upgraded["coordinationSpace"]["obsLabelsType"] == {"A": "existing", "B": "a", "C": "b"}
    assert upgraded["layout"][0]["coordinationScopes"] == {"dataset": "A", "obsLabelsType": ["B", "C"]}
    assert upgraded["layout"][1]["coordinationScopes"] == {"obsLabelsType": ["B", "C"]}
    assert upgraded["layout"][2]["coordinationScopes"] == {"dataset": "B"}


def test_upgrade_from_1_0_12_without_coordination_space_or_datasets():
    config = files_config("1.0.12", {"fileType": "anndata-cells.zarr", "options": {"factors": ["obs/a", "obs/b"]}})
    assert upgrade_config(config)["coordinationSpace"]["obsLabelsType"] == {"A": "a", "B": "b"}
    assert upgrade_config({"version": "1.0.12", "datasets": [], "layout": [view("status")]})["coordinationSpace"] == {}


def test_upgrade_from_1_0_14():
    config = {"version": "1.0.14", "datasets": [], "layout": [view("genes", props={"variablesLabelOverride": "gene"}), view("heatmap")]}
    with pytest.warns(UserWarning, match="'variablesLabelOverride' prop on the featureList view is deprecated"):
        upgraded = upgrade_config(config)
    assert [v["component"] for v in upgraded["layout"]] == ["featureList", "heatmap"]


@pytest.mark.parametrize("old_key", ["labels", "shapes"])
def test_upgrade_from_1_0_17(old_key):
    config = files_config("1.0.17", {"fileType": "spatialdata.zarr", "options": {old_key: {"path": "x"}}}, {"fileType": "spatialdata.zarr.zip"})
    files = upgrade_config(config)["datasets"][0]["files"]
    assert files[0]["options"] == {"obsSegmentations": {"path": "x"}}
    assert "options" not in files[1]


def test_upgrade_from_1_0_17_labels_and_shapes():
    config = files_config("1.0.17", {"fileType": "spatialdata.zarr", "options": {"labels": {}, "shapes": {}}})
    with pytest.raises(ValueError, match="Cannot upgrade both labels and shapes"):
        upgrade_config(config)


def test_upgrade_unknown_version():
    with pytest.raises(ValueError, match="not recognized: '2.0.0'"):
        upgrade_config({"version": "2.0.0"})
