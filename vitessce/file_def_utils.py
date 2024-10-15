from __future__ import annotations

from functools import partial
from typing import Optional

import numpy as np


def gen_obs_embedding_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None, dims: Optional[list[list[int]]] = None):
    if paths is not None:
        if "obsEmbedding" not in options:
            options["obsEmbedding"] = []
        if names is not None:
            for key, mapping in zip(paths, names):
                options["obsEmbedding"].append({
                    "path": key,
                    "dims": [0, 1],
                    "embeddingType": mapping
                })
        else:
            for mapping in paths:
                mapping_key = mapping.split('/')[-1]
                options["obsEmbedding"].append({
                    "path": mapping,
                    "dims": [0, 1],
                    "embeddingType": mapping_key
                })
    if dims is not None:
        if "obsEmbedding" not in options:
            options["obsEmbedding"] = []
        for dim_i, dim in enumerate(dims):
            options["obsEmbedding"][dim_i]['dims'] = dim
    return options


def gen_obs_sets_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None):
    if paths is not None:
        options["obsSets"] = []
        if names is not None:
            names = names
        else:
            names = []
            for obs in paths:
                obs_end_path = obs.split('/')[-1]
                names += [obs_end_path]
        for obs, name in zip(paths, names):
            options["obsSets"].append({
                "name": name,
                "path": obs
            })
    return options


def gen_sdata_obs_sets_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None, table_path: Optional[str] = None, region: Optional[str] = None):
    if paths is not None:
        options["obsSets"] = {"obsSets": []}
        if names is not None:
            names = names
        else:
            names = []
            for obs in paths:
                obs_end_path = obs.split('/')[-1]
                names += [obs_end_path]
        for obs, name in zip(paths, names):
            options["obsSets"]["obsSets"].append({
                "name": name,
                "path": obs
            })
        if table_path is not None:
            options["obsSets"]["tablePath"] = table_path
        if region is not None:
            options["obsSets"]["region"] = region
    return options


def gen_obs_feature_matrix_schema(options: dict, matrix_path: Optional[str] = None, var_filter_path: Optional[str] = None, init_var_filter_path: Optional[str] = None):
    if matrix_path is not None:
        options["obsFeatureMatrix"] = {
            "path": matrix_path
        }
        if var_filter_path is not None:
            options["obsFeatureMatrix"]["featureFilterPath"] = var_filter_path
        if init_var_filter_path is not None:
            options["obsFeatureMatrix"]["initialFeatureFilterPath"] = init_var_filter_path
    return options


def gen_obs_labels_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None):
    if paths is not None:
        if names is not None and len(paths) == len(names):
            # A name was provided for each path element, so use those values.
            names = names
        else:
            # Names were not provided for each path element,
            # so fall back to using the final part of each path for the names.
            names = [labels_path.split('/')[-1] for labels_path in paths]
        obs_labels = []
        for path, name in zip(paths, names):
            obs_labels.append({"path": path, "obsLabelsType": name})
        options["obsLabels"] = obs_labels
    return options


def gen_path_schema(key: str, path: Optional[str], options: dict):
    if path is not None:
        options[key] = {
            "path": path
        }
    return options


gen_obs_locations_schema = partial(gen_path_schema, "obsLocations")
gen_obs_segmentations_schema = partial(gen_path_schema, "obsSegmentations")
gen_obs_spots_schema = partial(gen_path_schema, "obsSpots")
gen_obs_points_schema = partial(gen_path_schema, "obsPoints")
gen_feature_labels_schema = partial(gen_path_schema, "featureLabels")


def gen_sdata_image_schema(options, path: str, coordinate_system: Optional[str] = None, affine_transformation: Optional[np.ndarray] = None) -> dict:
    if path is not None:
        options["image"] = {
            "path": path
        }
        if affine_transformation is not None:
            options["image"]['coordinateTransformations'] = affine_transformation
        if coordinate_system is not None:
            options["image"]['coordinateSystem'] = coordinate_system
    return options


def gen_sdata_labels_schema(options, path: str, table_path: str = "tables/table", coordinate_system: Optional[str] = None, affine_transformation: Optional[np.ndarray] = None) -> dict:
    if path is not None:
        options["labels"] = {
            "path": path
        }
        if table_path is not None:
            options["labels"]['tablePath'] = table_path
        if affine_transformation is not None:
            options["labels"]['coordinateTransformations'] = affine_transformation
        if coordinate_system is not None:
            options["labels"]['coordinateSystem'] = coordinate_system
    return options


def gen_sdata_obs_spots_schema(options: dict, shapes_path: str, table_path: str = "tables/table", region: Optional[str] = None, coordinate_system: Optional[str] = None) -> dict:
    if shapes_path is not None:
        options['obsSpots'] = {
            "path": shapes_path,
            "tablePath": table_path
        }
        if region is not None:
            options['obsSpots']['region'] = region
        if coordinate_system is not None:
            options['obsSpots']['coordinateSystem'] = coordinate_system
    return options


def gen_sdata_obs_feature_matrix_schema(options: dict, matrix_path: Optional[str] = None, var_filter_path: Optional[str] = None, init_var_filter_path: Optional[str] = None, region: Optional[str] = None):
    if matrix_path is not None:
        options["obsFeatureMatrix"] = {
            "path": matrix_path
        }
        if region is not None:
            options['obsFeatureMatrix']['region'] = region
        if var_filter_path is not None:
            options["obsFeatureMatrix"]["featureFilterPath"] = var_filter_path
        if init_var_filter_path is not None:
            options["obsFeatureMatrix"]["initialFeatureFilterPath"] = init_var_filter_path
    return options
