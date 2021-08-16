
## [0.1.0a13] (In Progress)

### Added
- Added a `VitessceConfig.web_app()` method for launching the `vitessce.io` web application based on the current config object, and serving the data associated with the config via the Starlette server.

### Changed
- Bump version of config to 1.0.3.  `disableChannelsIfRgbDetected` prop now needs to be used on the `LayerController` component to hide the channel controller for detected RGB images.  Fix notebooks for this.

## [0.1.0a12]

### Added
- Added a `VitessceConfig.widget()` method for convenience when creating the widget (no more need to import a separate `VitessceWidget` class).
- Added a `VitessceConfig.link_views()` method for convenience when adding coordinations and setting default coordination scope values.
- Added `pandas` dependency.

### Changed
- Removed old notebooks.


## [0.1.0a3](https://pypi.org/project/vitessce/0.1.0a3/)
