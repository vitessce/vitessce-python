# Example Notebooks

The notebooks contained in this directory demonstrate the Vitessce Python package functionality.

## Setup

Some of the example notebooks rely on external single-cell data analysis packages. An environment containing these additional packages can be installed with `conda` or `pip`.

```sh
cd docs/notebooks
uv sync --extra dev --extra docs --extra all
```

## Run

```sh
uv run jupyter lab
```

## Troubleshooting

If you previously had the `vitessce` Python package v1 installed, you may need to uninstall the previous lab extension due to conflicts between the v1 and v2 widget JS code.

```sh
jupyter labextension uninstall vitessce-jupyter
```
