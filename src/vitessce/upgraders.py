"""
Upgrade view configs from older schema versions to the latest one.

This is a port of ``previous-config-upgraders.js`` from the ``@vitessce/schemas``
NPM package (version 4.1.3), with one function per JS ``upgradeFrom*`` function.
JS upgraders that only bump the version have no Python function.
The ``_upgrade_from_*`` functions mutate and return their argument;
``upgrade_config`` deep-copies its input first.
"""
import copy
import uuid
import warnings

from .config import _get_next_scope

VERSIONS = ["0.1.0"] + [f"1.0.{i}" for i in range(19)]


def _js_truthy(value):
    # Unlike Python, JS treats [] and {} as truthy.
    return value not in (None, False, 0, "")


def _replace_view_prop(prefix, view, coordination_space):
    values = {"Zoom": view["zoom"], "TargetX": view["target"][0], "TargetY": view["target"][1]}
    new_scopes = {}
    for suffix, value in values.items():
        c_type = f"{prefix}{suffix}"
        scope = _get_next_scope(list(coordination_space[c_type].keys()))
        coordination_space[c_type][scope] = value
        new_scopes[c_type] = scope
    return new_scopes


def _upgrade_from_0_1_0(config, dataset_uid=None):
    coordination_space = {c_type: {} for c_type in [
        "embeddingType", "embeddingZoom", "embeddingTargetX", "embeddingTargetY",
        "spatialZoom", "spatialTargetX", "spatialTargetY",
    ]}
    layout = []
    for component_def in config["staticLayout"]:
        new_def = {**component_def, "coordinationScopes": {}}
        props = component_def.get("props") or {}
        if component_def["component"] == "scatterplot":
            mapping = props.get("mapping")
            if mapping and isinstance(mapping, str):
                coordination_space["embeddingType"][mapping] = mapping
                new_def["coordinationScopes"]["embeddingType"] = mapping
            if _js_truthy(props.get("view")):
                new_def["coordinationScopes"].update(_replace_view_prop("embedding", props["view"], coordination_space))
        if component_def["component"] == "spatial" and _js_truthy(props.get("view")):
            new_def["coordinationScopes"].update(_replace_view_prop("spatial", props["view"], coordination_space))
        layout.append(new_def)

    lc_def = next((c for c in layout if c["component"] == "layerController"), None)
    spatial_def = next((c for c in layout if c["component"] == "spatial"), None)
    if lc_def and spatial_def:
        lc_def["coordinationScopes"] = dict(spatial_def["coordinationScopes"])

    dataset_uid = dataset_uid or str(uuid.uuid4())
    return {
        "version": "1.0.0",
        **{key: config[key] for key in ("name", "description", "public") if key in config},
        "datasets": [{
            "uid": dataset_uid,
            "name": dataset_uid,
            "files": [{"fileType": layer["fileType"], "url": layer["url"]} for layer in config["layers"]],
        }],
        "initStrategy": "auto",
        "coordinationSpace": coordination_space,
        "layout": layout,
    }


_SPATIAL_LAYER_TYPES = {
    "raster": "spatialRasterLayers",
    "cells": "spatialCellsLayer",
    "molecules": "spatialMoleculesLayer",
    "neighborhoods": "spatialNeighborhoodsLayer",
}


def _upgrade_from_1_0_0(config):
    coordination_space = config.get("coordinationSpace")
    if coordination_space and "spatialLayers" in coordination_space:
        for layer_type, c_type in _SPATIAL_LAYER_TYPES.items():
            coordination_space[c_type] = {}
            for scope, layers in coordination_space["spatialLayers"].items():
                typed_layers = [
                    {k: v for k, v in layer.items() if k != "type"}
                    for layer in (layers if isinstance(layers, list) else [])
                    if layer.get("type") == layer_type
                ]
                if not typed_layers:
                    coordination_space[c_type][scope] = None
                else:
                    coordination_space[c_type][scope] = typed_layers if layer_type == "raster" else typed_layers[0]
        del coordination_space["spatialLayers"]

    for component in config["layout"]:
        scopes = component.get("coordinationScopes")
        if scopes and _js_truthy(scopes.get("spatialLayers")):
            for layer_type, c_type in _SPATIAL_LAYER_TYPES.items():
                if component["component"] in ("spatial", "layerController") or (component["component"] == "description" and layer_type == "raster"):
                    scopes[c_type] = scopes["spatialLayers"]
            del scopes["spatialLayers"]
    return config


def _set_layer_controller_prop(config, prop):
    for component in config["layout"]:
        if component["component"] == "layerController":
            component["props"] = {**(component.get("props") or {}), prop: True}


def _upgrade_from_1_0_1(config):
    # Match the previous lack of 3D auto-detection behavior.
    _set_layer_controller_prop(config, "globalDisable3d")
    # Enforce bitmask or raster as the spatial raster layer type.
    raster_layers = (config.get("coordinationSpace") or {}).get("spatialRasterLayers") or {}
    for layers in raster_layers.values():
        for layer in layers or []:
            if layer.get("type") not in ("bitmask", "raster"):
                layer["type"] = "raster"
    return config


def _upgrade_from_1_0_2(config):
    _set_layer_controller_prop(config, "disableChannelsIfRgbDetected")
    return config


_COORDINATION_TYPE_ANALOGIES_1_0_11 = {
    # Spatial layer types
    "spatialRasterLayers": "spatialImageLayer",
    "spatialCellsLayer": "spatialSegmentationLayer",
    "spatialMoleculesLayer": "spatialPointLayer",
    "spatialNeighborhoodsLayer": "spatialNeighborhoodLayer",
    # Other types
    "cellFilter": "obsFilter",
    "cellHighlight": "obsHighlight",
    "cellSelection": "obsSelection",
    "cellSetSelection": "obsSetSelection",
    "cellSetHighlight": "obsSetHighlight",
    "cellSetColor": "obsSetColor",
    "geneFilter": "featureFilter",
    "geneHighlight": "featureHighlight",
    "geneSelection": "featureSelection",
    "geneExpressionColormap": "featureValueColormap",
    "geneExpressionColormapRange": "featureValueColormapRange",
    "cellColorEncoding": "obsColorEncoding",
    "additionalCellSets": "additionalObsSets",
    "embeddingCellSetPolygonsVisible": "embeddingObsSetPolygonsVisible",
    "embeddingCellSetLabelsVisible": "embeddingObsSetLabelsVisible",
    "embeddingCellSetLabelSize": "embeddingObsSetLabelSize",
    "embeddingCellRadius": "embeddingObsRadius",
    "embeddingCellRadiusMode": "embeddingObsRadiusMode",
    "embeddingCellOpacity": "embeddingObsOpacity",
    "embeddingCellOpacityMode": "embeddingObsOpacityMode",
}


def _rename_truthy_keys(obj, analogies):
    for old_key, new_key in analogies.items():
        if _js_truthy(obj.get(old_key)):
            obj[new_key] = obj.pop(old_key)


def _upgrade_from_1_0_10(config):
    # Changes to spatial layer coordination type names; Cell -> Obs, Gene -> Feature.
    config["coordinationSpace"] = config.get("coordinationSpace") or {}
    _rename_truthy_keys(config["coordinationSpace"], _COORDINATION_TYPE_ANALOGIES_1_0_11)
    for component in config["layout"]:
        component["coordinationScopes"] = component.get("coordinationScopes") or {}
        _rename_truthy_keys(component["coordinationScopes"], _COORDINATION_TYPE_ANALOGIES_1_0_11)
    return config


def _upgrade_from_1_0_11(config):
    # Pass the embedding types to cells.json files, since the fileType-to-dataType mapping replaced files[].type.
    coordination_space = config.get("coordinationSpace") or {}
    if _js_truthy(coordination_space.get("embeddingType")):
        embedding_types = list(coordination_space["embeddingType"].values())
        for dataset in config["datasets"]:
            for file_def in dataset["files"]:
                if file_def["fileType"] == "cells.json":
                    file_def["options"] = {"embeddingTypes": list(embedding_types)}
    return config


def _upgrade_from_1_0_12(config):
    # Set up obsLabelsType coordination scopes for anndata-cells.zarr options.factors.
    datasets = config["datasets"]
    coordination_space = config.get("coordinationSpace") or {}
    dataset_obs_labels_type_scopes = {}
    for dataset in datasets:
        for file_def in dataset["files"]:
            options = file_def.get("options")
            if file_def["fileType"] == "anndata-cells.zarr" and isinstance(options, dict) and isinstance(options.get("factors"), list):
                scopes = []
                for obs_labels_type in options["factors"]:
                    # JS reads the pre-upgrade space here, which repeats scope "A" when coordinationSpace is missing.
                    next_scope = _get_next_scope(list(coordination_space.get("obsLabelsType", {}).keys()))
                    # Remove the obs/ prefix.
                    coordination_space.setdefault("obsLabelsType", {})[next_scope] = obs_labels_type.split("/")[-1]
                    scopes.append(next_scope)
                dataset_obs_labels_type_scopes[dataset["uid"]] = scopes

    for view in config["layout"]:
        dataset_scope = (view.get("coordinationScopes") or {}).get("dataset")
        if dataset_scope and isinstance(dataset_scope, str):
            view_dataset_uid = coordination_space.get("dataset", {}).get(dataset_scope)
        else:
            view_dataset_uid = datasets[0]["uid"] if datasets else None
        if isinstance(view_dataset_uid, str) and view_dataset_uid in dataset_obs_labels_type_scopes:
            view["coordinationScopes"] = {
                **(view.get("coordinationScopes") or {}),
                "obsLabelsType": list(dataset_obs_labels_type_scopes[view_dataset_uid]),
            }
    config["coordinationSpace"] = coordination_space
    return config


_VIEW_TYPE_ANALOGIES_1_0_15 = {
    "genes": "featureList",
    "cellSets": "obsSets",
    "cellSetSizes": "obsSetSizes",
    "cellSetExpression": "obsSetFeatureValueDistribution",
    "expressionHistogram": "featureValueHistogram",
}

_PROP_ANALOGIES_1_0_15 = {
    "variablesLabelOverride": "featureType",
    "observationsLabelOverride": "obsType",
}


def _upgrade_from_1_0_14(config):
    # View type renames, and deprecation of the label override props.
    for view in config["layout"]:
        view["component"] = _VIEW_TYPE_ANALOGIES_1_0_15.get(view["component"], view["component"])
        for old_prop, new_type in _PROP_ANALOGIES_1_0_15.items():
            if _js_truthy((view.get("props") or {}).get(old_prop)):
                warnings.warn(f"The '{old_prop}' prop on the {view['component']} view is deprecated. Please use the '{new_type}' coordination type instead.")
    return config


def _upgrade_from_1_0_15(config):
    # Move implicit per-dataset coordination scopes into the explicit coordinationScopesBy property.
    for view in config["layout"]:
        scopes = view.get("coordinationScopes") or {}
        view["coordinationScopes"] = {c_type: scope for c_type, scope in scopes.items() if isinstance(scope, (list, str))}
        if isinstance(scopes.get("dataset"), list):
            view["coordinationScopesBy"] = {
                "dataset": {c_type: scope for c_type, scope in scopes.items() if isinstance(scope, dict)},
            }
    return config


def _upgrade_from_1_0_16(config):
    # Convert obsSets.anndata.zarr and obsFeatureColumns.anndata.zarr options from array into object property.
    option_keys = {"obsSets.anndata.zarr": "obsSets", "obsFeatureColumns.anndata.zarr": "obsFeatureColumns"}
    for dataset in config["datasets"]:
        for file_def in dataset["files"]:
            if file_def["fileType"] in option_keys:
                file_def["options"] = {option_keys[file_def["fileType"]]: file_def.get("options")}
    return config


def _upgrade_from_1_0_17(config):
    # Rename spatialdata.zarr labels/shapes options to obsSegmentations.
    for dataset in config["datasets"]:
        for file_def in dataset["files"]:
            if file_def["fileType"] in ("spatialdata.zarr", "spatialdata.zarr.zip"):
                options = file_def.get("options") or {}
                if _js_truthy(options.get("labels")) and _js_truthy(options.get("shapes")):
                    raise ValueError("Cannot upgrade both labels and shapes at the same time. For now, a workaround is to include two separate spatialdata.zarr file definitions, one for labels and one for shapes.")
                for old_key in ("labels", "shapes"):
                    if _js_truthy(options.get(old_key)):
                        options["obsSegmentations"] = options.pop(old_key)
    return config


UPGRADERS = {
    "0.1.0": _upgrade_from_0_1_0,
    "1.0.0": _upgrade_from_1_0_0,
    "1.0.1": _upgrade_from_1_0_1,
    "1.0.2": _upgrade_from_1_0_2,
    "1.0.10": _upgrade_from_1_0_10,
    "1.0.11": _upgrade_from_1_0_11,
    "1.0.12": _upgrade_from_1_0_12,
    "1.0.14": _upgrade_from_1_0_14,
    "1.0.15": _upgrade_from_1_0_15,
    "1.0.16": _upgrade_from_1_0_16,
    "1.0.17": _upgrade_from_1_0_17,
}


def upgrade_config(config):
    """
    Upgrade a view config to the latest schema version, as the Vitessce JS app does when it loads a config.
    Check the config with ``validate_config`` first, since the upgrade steps assume a valid config.

    :param dict config: A view config.
    :returns: An upgraded copy of the config.
    :rtype: dict
    :raises ValueError: If the version is not recognized, or if an upgrade step cannot convert the config.

    .. code-block:: python

        from vitessce import validate_config, upgrade_config

        if not validate_config(config):
            config = upgrade_config(config)
    """
    version = config.get("version")
    if version not in VERSIONS:
        raise ValueError(f"Config version was not recognized: {version!r}")
    config = copy.deepcopy(config)
    i = VERSIONS.index(version)
    for prev_version, next_version in zip(VERSIONS[i:], VERSIONS[i + 1:]):
        if prev_version in UPGRADERS:
            config = UPGRADERS[prev_version](config)
        config["version"] = next_version
    return config
