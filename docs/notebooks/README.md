# Example Notebooks

The notebooks contained in this directory demonstrate the Vitessce Python package functionality.

## Setup

Some of the example notebooks rely on external single-cell data analysis packages. An environment containing these additional packages can be installed with `conda` or `pip`.

```sh
cd docs/notebooks
conda env create -f environment.yml
conda activate vitessce-python-notebooks
pip install -e "../..[dev]"
jupyter labextension install ../../js
```

## Run

```sh
jupyter lab
```
