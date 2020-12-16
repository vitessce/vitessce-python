# Example Notebooks

The notebooks contained in this directory demonstrate the Vitessce Jupyter widget functionality.

## Setup

Some of the example notebooks rely on external single-cell data analysis packages. An environment containing these additional packages can be installed with `conda` or `pip`.

```sh
conda env create -f environment.yml
conda activate vitessce-jupyter-examples
pip install -e ..
jupyter nbextension install --py --symlink --sys-prefix vitessce
jupyter nbextension enable --py --sys-prefix vitessce
```

## Run

```sh
jupyter notebook
```