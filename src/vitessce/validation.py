"""
Validate view configs against the JSON schemas from the ``@vitessce/json-schema`` NPM package.

The ``json_schema/`` directory holds the ``dist/`` files of that package at ``JSON_SCHEMA_VERSION``.
To refresh them, for example when bumping ``DEFAULT_JS_PACKAGE_VERSION`` in ``widget.py``, run:

    curl -sL https://registry.npmjs.org/@vitessce/json-schema/-/json-schema-<version>.tgz | tar xz -C /tmp
    cp /tmp/package/dist/*.json src/vitessce/json_schema/
"""
import json
from functools import lru_cache
from importlib.resources import files

from jsonschema import Draft7Validator

from .upgraders import VERSIONS, upgrade_config

JSON_SCHEMA_VERSION = "4.1.3"
LATEST_VERSION = VERSIONS[-1]


@lru_cache
def _validator(name):
    return Draft7Validator(json.loads(files("vitessce").joinpath("json_schema", name).read_text()))


def _is_never_match(error):
    # A const on a direct child is a discriminator (such as fileType) that selects a different anyOf branch.
    # A root-level `not: {}` is how zod-to-json-schema encodes the "undefined" branch of an optional value.
    return (error.validator == "const" and len(error.relative_path) == 1) or (error.validator == "not" and not error.relative_path)


def _is_type_mismatch(error):
    return error.validator == "type" and not error.relative_path


def _leaf_errors(error):
    # Descend into an anyOf error, following the one branch that was meant for the instance.
    if not error.context:
        return [{"path": error.json_path, "message": error.message}]
    branches = {}
    for sub_error in error.context:
        branches.setdefault(sub_error.schema_path[0], []).append(sub_error)
    kept = [branch for branch in branches.values() if not any(_is_never_match(sub_error) for sub_error in branch)]
    if len(kept) > 1:
        # Prefer branches for the instance's JSON type, such as the array branch of a nullable array.
        kept = [branch for branch in kept if not any(_is_type_mismatch(sub_error) for sub_error in branch)] or kept
    if len(kept) == 1:
        return [leaf for sub_error in kept[0] for leaf in _leaf_errors(sub_error)]
    consts = [sub_error for sub_error in error.context if sub_error.validator == "const"]
    if not kept and consts:
        return [{"path": consts[0].json_path, "message": f"{consts[0].instance!r} is not a known {consts[0].relative_path[0]}"}]
    return [{"path": error.json_path, "message": error.message}]


def _errors(schema_name, config):
    return [leaf for error in _validator(schema_name).iter_errors(config) for leaf in _leaf_errors(error)]


def validate_config(config, base_plugins=True):
    """
    Validate a view config, as the Vitessce JS app does when it loads a config.
    A config with an older version is checked against the schema for that version,
    upgraded with ``upgrade_config``, then checked against the latest schema.

    :param dict config: A view config, such as the output of ``VitessceConfig.to_dict()``.
    :param bool base_plugins: Whether to check file options, view components, and coordination values
        against the plugins that ship with Vitessce. Set to False for configs that use plugin view types,
        file types, or coordination types. By default, True.
    :returns: The errors. Each is a dict with a JSONPath ``"path"`` and a ``"message"``.
        Paths of errors found after the upgrade refer to the upgraded config.
        The list is empty when the config is valid.
    :rtype: list[dict]

    .. code-block:: python

        from vitessce import validate_config

        errors = validate_config(vc.to_dict(base_url="http://localhost:8000"))
        # [{'path': '$.layout[0].component', 'message': "'scatterplott' is not one of [...]"}]
    """
    if not isinstance(config, dict):
        return [{"path": "$", "message": "config must be an object"}]
    version = config.get("version")
    if version not in VERSIONS:
        return [{"path": "$.version", "message": f"Config version was not recognized: {version!r}"}]
    if version != LATEST_VERSION:
        errors = _errors(f"config-{version}.schema.json", config)
        if errors:
            return errors
        try:
            config = upgrade_config(config)
        except ValueError as e:
            return [{"path": "$", "message": str(e)}]
    schema_name = f"config-{LATEST_VERSION}.base-plugins.schema.json" if base_plugins else f"config-{LATEST_VERSION}.schema.json"
    return _errors(schema_name, config)
