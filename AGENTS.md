# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## Commands

```sh
uv sync --extra dev --extra docs --extra all   # dev install (add --extra notebook for marimo/spatialdata notebooks)
uv run flake8                                  # lint (config in setup.cfg; E501 line length ignored)
uv run pytest                                  # all tests (testpaths = tests/)
uv run pytest tests/test_config.py::test_config_creation   # single test
uv run autopep8 --in-place --aggressive path/to/file.py    # format
uv run make html                               # Sphinx docs -> docs-dist/html (needs pandoc)
```

CI coverage gate (`.github/workflows/test.yml`):

```sh
uv run coverage run --rcfile .coveragerc_omit --module pytest && uv run coverage report --show-missing --fail-under 100
```

Every module under `src/` not listed in `.coveragerc_omit` must have 100% coverage. A new module needs full test coverage or an entry in that omit list.

Widget render tests (`tests-widget/`, Playwright + pnpm) run against HTML exported from notebooks:

```sh
uv run jupyter nbconvert --to=html --execute docs/notebooks/widget_from_dict.ipynb
uv run marimo export html docs/notebooks/marimo.py -o docs/notebooks/marimo.html
cd tests-widget && pnpm install && pnpm exec playwright test
```

Demo pipelines: `cd demos && uv run snakemake -j 1 --dry-run` (needs `--extra demos`).

## Architecture

The package builds Vitessce view configs (JSON) in Python and renders them with the Vitessce JS app, either in a notebook widget or on vitessce.io.

- **`config.py`**: `VitessceConfig` is the builder. It handles datasets, views, coordination scopes, and layout (`hconcat`/`vconcat`, `|` and `/` operators). `VitessceChainableConfig` is the fluent variant that `to_python()` emits. `to_dict(base_url=...)` produces the final JSON. Data URLs resolve only at that point, so one config can target a localhost server, S3 (`export_to_s3`), or static files (`export_to_files`).
- **`constants.py`**: `ViewType`, `CoordinationType`, `DataType`, and `FileType` enums. They mirror the string constants in the JS Vitessce package. `norm_enum` accepts either the enum or a raw string, so plugin-defined types also work.
- **`wrappers.py`**: `AbstractWrapper` subclasses (`AnnDataWrapper`, `SpatialDataWrapper`, OME-TIFF/OME-Zarr wrappers, etc.) adapt local or remote data into Vitessce file definitions. Each wrapper's `convert_and_save(dataset_uid, obj_i)` fills three things:
  - `file_def_creators`: functions of `base_url` that return file-definition dicts.
  - `routes`: Starlette routes for local files served over HTTP.
  - `stores`: zarr stores served over the widget comm channel instead of HTTP.

  Wrappers accept exactly one of `*_path`, `*_url`, `*_store`, or `*_artifact` (lamindb). The `raise_error_if_*` helpers enforce this. Each wrapper sets `self._repr = make_repr(locals())` in `__init__` so `VitessceConfig.to_python()` can regenerate constructor calls.
- **`widget.py`**: `VitessceWidget` is an `anywidget.AnyWidget`. The front end is the inline `ESM` string, which loads the `vitessce` npm package from unpkg at `js_package_version`. Python serves data to the front end in two ways:
  1. A uvicorn `BackgroundServer` on port 8000 and up serves wrapper `routes` (`serve_routes`).
  2. Zarr `stores` are read through custom comm messages (`_zarr_get`, `_zarr_get_range`, `_zarr_get_multi`, dispatched in `_dispatch_command`). These handlers must stay synchronous; `sync_store.py` wraps async zarr v3 stores.

  `VitesscePlugin` adds `plugin_esm`, RPC `commands` (via `_plugin_command`), and an `on_config_change` hook. Example plugins live in `widget_plugins/`. `ipython_display` is a non-ipywidgets HTML fallback that reuses the same ESM.
- **`data_utils/`**: helpers that convert data into Vitessce-friendly formats (AnnData optimization, OME-TIFF/Zarr writing, multivec, SpatialData points Morton z-order).
- **`config_converter.py`**: converts UCSC Cell Browser projects to AnnData-Zarr plus a config.

Optional dependencies: the base install has only scipy, pandas, numpy, zarr, black, and jsonschema. `widget`, `wrappers`, and `export` need `[all]`. `__init__.py` wraps those imports in `try/except ModuleNotFoundError`. Several heavy imports (starlette, widget) are deferred into function bodies so `import vitessce` works without extras. Keep that pattern in modules the base install imports.

## Releasing

- **Python package:** bump `version` in `pyproject.toml`. A merge to `main` triggers `deploy.yml`, which publishes to PyPI and deploys the docs to GitHub Pages.
- **Bundled Vitessce JS version:** update `DEFAULT_JS_PACKAGE_VERSION` in `src/vitessce/widget.py`. Also refresh the vendored schemas in `src/vitessce/json_schema/` to the same `@vitessce/json-schema` version, using the command in the `validation.py` docstring, and update `JSON_SCHEMA_VERSION`. `tests/test_validation.py` fails if the two versions differ.
- **New view config schema version:** port any new `upgradeFrom*` function from `@vitessce/schemas` `previous-config-upgraders.js` into `src/vitessce/upgraders.py`. Then extend `VERSIONS` (which also sets `LATEST_VERSION`) and `UPGRADERS`.
