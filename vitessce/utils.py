from functools import partial
from typing import Optional

import numpy as np


def get_next_scope_numeric(prev_scopes):
    next_scope_int = 0
    next_scope_str = None

    while True:
        next_scope_str = str(next_scope_int)
        if next_scope_str not in prev_scopes:
            break
        next_scope_int += 1
    return next_scope_str


def create_prefixed_get_next_scope_numeric(prefix):

    def inner_get_next_scope(prev_scopes):
        next_scope_int = 0
        next_scope_str = None

        while True:
            next_scope_str = f"{prefix}{next_scope_int}"
            if next_scope_str not in prev_scopes:
                break
            next_scope_int += 1
        return next_scope_str

    return inner_get_next_scope


def get_initial_coordination_scope_prefix(dataset_uid, data_type):
    return f"init_{dataset_uid}_{data_type}_"


def get_initial_coordination_scope_name(dataset_uid, data_type, i=None):
    prefix = get_initial_coordination_scope_prefix(dataset_uid, data_type)
    return f"{prefix}{0 if i is None else i}"



def gen_obs_embedding_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None, dims: Optional[list[list[int]]] = None):
    if paths is not None:
        if names is not None:
            for key, mapping in zip(paths, names):
                options["obsEmbedding"].append({
                    "path": mapping,
                    "dims": [0, 1],
                    "embeddingType": key
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
        for dim_i, dim in enumerate(dims):
            options["obsEmbedding"][dim_i]['dims'] = dim
    return options

def gen_obs_sets_schema(options: dict, paths: Optional[list[str]] = None, names: Optional[list[str]] = None):
    if paths is not None:
        options["obsSets"] = []
        if names is not None:
            names = names
        else:
            names = [obs.split('/')[-1] for obs in paths]
        for obs, name in zip(paths, names):
            options["obsSets"].append({
                "name": name,
                "path": obs
            })
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

gen_obs_locations_schema = partial(gen_path_schema,  "obsLocations")
gen_obs_segmentations_schema = partial(gen_path_schema,  "obsSegmentations")
gen_obs_spots_schema = partial(gen_path_schema,  "obsSpots")
gen_obs_points_schema = partial(gen_path_schema,  "obsPoints")
gen_feature_labels_schema = partial(gen_path_schema, "featureLabels")

def gen_image_schema(options, path: str, affine_transformation: Optional[np.ndarray] = None) -> dict:
    if path is not None:
        options["image"] = {
            "path": path
        }
        if affine_transformation is not None:
            options['coordinateTransformations'] = affine_transformation
    return options

def gen_obs_spots_schema(options: dict, shapes_path: Optional[str] = None) -> dict:
    if shapes_path is not None:
        options['obsSpots'] = {
            "path": shapes_path,
            "tablePath": "table/table",
            "region": "region"
        }
    return options