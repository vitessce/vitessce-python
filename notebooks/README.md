# Example Notebooks

The notebooks contained here demonstrate the Vitessce Jupyter widget with different use cases.

- [Minimal Example](./example-minimal.ipynb) - Initialization of the Vitessce widget with a pre-defined view config.
- [Creation of a Custom Vitessce Config](./example-config-creation.ipynb) - Creation of a Vitessce config for displaying a local single-cell dataset.
- [Creation of a Custom Vitessce Config with Imaging Data](./example-config-creation-with-images.ipynb) - Creation of a Vitessce config for displaying a local single-cell dataset with associated microscopy images.
- [Creation of a Custom Vitessce Config with Coordinated Views](./example-config-creation-with-coordination.ipynb) - Creation of a Vitessce config that demonstrates the coordination model.
- [Cell Selection](./example-cell-selection.ipynb) - Use the Vitessce widget to make a selection of cells, then access the list of selected cell IDs in Python.
- [Differential Expression](./example-differential-expression.ipynb) - Use the Vitessce widget to make two cell selections, then compare the two selections by performing a differential expression analysis in Python.
- [Spatial Differential Expression](./example-differential-expression-spatial.ipynb) - Use the Vitessce widget to make two cell selections, then compare the two selections by performing a differential expression analysis in Python.


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