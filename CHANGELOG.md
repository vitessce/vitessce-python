## 3.0.7 - 2023-08-11
### Changed
- Added support for multiple coordination scopes of the same type using the `allow_multiple_scopes_per_type` named arg
- Added new `CoordinationType`s and `FileType`s
- Added support for `obs_labels_names` and `obs_labels_paths` in AnndataWrapper

## 1.0.10

### Changed
- Pin `ipywidgets` and `jupyterlab-widgets` dependency versions as a short-term fix. In the longer term, we will need to follow the [migration guide](https://ipywidgets.readthedocs.io/en/stable/migration_guides.html#migrating-from-7-x-to-8-0).

## 1.0.5 - in progress

### Added
- `VitessceConfig.to_python` method to auto-generate Python code representing a view config, intended to be used in auto-generated Jupyter notebooks

### Changed
- Simplify CI by removing bash CLI flags.
- Use `nvm` and `.nvmrc`.
- Update package-lock.json; Had gotten out of sync.
- Bump vitessce JS version to `1.1.17`
- Updated CHANGELOG
- Fix typo in constants; Should error if the same mistake is made in the future.


## [1.0.4](https://pypi.org/project/vitessce/1.0.4/) - 2021-09-28

### Changed
- Updated `setup.py` and README installation/development instructions based on latest [ipywidget cookiecutter template](https://github.com/jupyter-widgets/widget-cookiecutter/tree/9694718)


## [0.1.0a13](https://pypi.org/project/vitessce/0.1.0a13/) - 2021-09-11

### Added
- Added a `VitessceConfig.web_app()` method for launching the `vitessce.io` web application based on the current config object, and serving the data associated with the config via the Starlette server.

### Changed
- Bump version of config to 1.0.3.  `disableChannelsIfRgbDetected` prop now needs to be used on the `LayerController` component to hide the channel controller for detected RGB images.  Fix notebooks for this.
- Bump version of config to 1.0.4.
- Added a check for `divRef.current` being null in the cleanup function in `js/lib/widget.js`.


## [0.1.0a12](https://pypi.org/project/vitessce/0.1.0a12/) - 2021-07-20

### Added
- Added a `VitessceConfig.widget()` method for convenience when creating the widget (no more need to import a separate `VitessceWidget` class).
- Added a `VitessceConfig.link_views()` method for convenience when adding coordinations and setting default coordination scope values.
- Added `pandas` dependency.

### Changed
- Removed old notebooks.


## [0.1.0a0](https://pypi.org/project/vitessce/0.1.0a3/) - 2020-11-19

- First version pushed to PyPI
