# vitessce-python

[![PyPI](https://img.shields.io/pypi/v/vitessce)](https://pypi.org/project/vitessce)
[![docs](https://img.shields.io/badge/docs-📖-57B4E9.svg)](https://vitessce.github.io/vitessce-python/)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vitessce/vitessce-python/main?filepath=docs/notebooks/widget_pbmc.ipynb)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vitessce/vitessce-python/blob/main/docs/notebooks/widget_on_colab.ipynb)


Python API and Jupyter widget facilitating interactive visualization of spatial single-cell data with [Vitessce](https://github.com/vitessce/vitessce).


## Installation

To install with pip:

    $ pip install vitessce[all]

## Getting started

Explore our [example notebooks](./docs/notebooks/).
These contain demos of different use cases and integrations with single-cell data analysis packages.


## Development

For a development installation (requires NodeJS and NPM),

    $ git clone https://github.com/vitessce/vitessce-python.git
    $ cd vitessce-python
    $ conda env create -f environment.yml
    $ conda activate vitessce-python-dev
    $ pip install -e ".[dev]"

### Conda environments

In this repository, there are multiple conda environments for different purposes:

- `vitessce-python-dev` (defined in [environment.yml](./environment.yml)) is used for the development of the `vitessce` package itself
- `vitessce-python-notebooks` (defined in [docs/notebooks/environment.yml](./docs/notebooks/environment.yml)) is used for running the example notebooks in the `docs/notebooks/` directory (see [`docs/notebooks/README.md`](./docs/notebooks#readme) for more information)
- `vitessce-python-binder` (defined in [binder/environment.yml](./binder/environment.yml)) is the environment used by Binder upon opening notebooks from this repository

## Linting and testing

```sh
flake8
pytest
```

### Formatting

One file:

```sh
autopep8 --in-place --aggressive ./path/to/file.py
```

All `.py` files:

```sh
find . -name '*.py' -exec autopep8 --in-place --aggressive '{}' \;
```


## Documentation

```sh
make html
```


## Deployment

To deploy a new version, increment the version of the Python package in [`setup.py`](./setup.py).

Then, when you push or merge the code with the incremented versions to `main`, the GitHub Action `deploy.yml` workflow will build and push the package to PyPI.


## Resources

- [ipywidget docs: Building a Custom Widget](https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20Custom.html)
- [ipywidget docs: Low Level Widget Tutorial](https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Low%20Level.html)
- [ipywidget example: ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet)
- [ipywidget example (with React): ipymaterialui](https://github.com/maartenbreddels/ipymaterialui)
- [ipywidget example (with React): higlass-python](https://github.com/higlass/higlass-python)
- [ipywidget cookiecutter](https://github.com/jupyter-widgets/widget-cookiecutter)
- [Sphinx: Getting Started](https://www.sphinx-doc.org/en/master/usage/quickstart.html)
- [Read the Docs Sphinx Theme](https://github.com/readthedocs/sphinx_rtd_theme)
- [jupyter server proxy](https://jupyter-server-proxy.readthedocs.io/en/latest/arbitrary-ports-hosts.html)

## Getting/Offering Help

If you have a specific bug or feature request, please feel free to open an [issue](https://github.com/vitessce/vitessce-python/issues/new).  Otherwise our [discussions](https://github.com/vitessce/vitessce-python/discussions) section is a great place to get help or offer it.  If you aren't sure if something is a bug or not, don't have all the reproduction steps, or just have a general question, feel free to open a discussion post.
