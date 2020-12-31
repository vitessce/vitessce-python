# vitessce-python

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vitessce/vitessce-python/master?filepath=docs/notebooks/widget_pbmc.ipynb)

Python API and Jupyter widget facilitating interactive visualization of spatial single-cell data with [Vitessce](https://github.com/vitessce/vitessce).


## Installation

To install with pip:

    $ pip install vitessce
    $ jupyter nbextension install --py vitessce
    $ jupyter nbextension enable --py vitessce

To enable the widget for Jupyter Lab run the following additional lines:

    $ jupyter labextension install @jupyter-widgets/jupyterlab-manager
    $ jupyter labextension install vitessce-jupyter


## Getting started

Explore our [example notebooks](./docs/notebooks/).
These contain demos of different use cases and integrations with single-cell data analysis packages.


## Development

For a development installation (requires npm),

    $ git clone https://github.com/vitessce/vitessce-python.git
    $ cd vitessce-python
    $ conda env create -f environment.yml
    $ conda activate vitessce-jupyter-dev
    $ pip install -e .
    $ jupyter nbextension install --py --symlink --sys-prefix vitessce
    $ jupyter nbextension enable --py --sys-prefix vitessce
    $ cd js
    $ jupyter labextension install @jupyter-widgets/jupyterlab-manager
    $ jupyter labextension install

When actively developing your extension, build Jupyter Lab with the command:

    $ jupyter lab --watch

This takes a minute or so to get started, but then automatically rebuilds JupyterLab when your javascript changes.

Note on first `jupyter lab --watch`, you may need to touch a file to get Jupyter Lab to open.


## Testing

```sh
cd tests
python -m unittest
```


## Documentation

```sh
make html
```


## Deployment

To deploy a new version, increment the version of the Python package in [`vitessce/_version.py`](./vitessce/_version.py) and the JS package in [`js/package.json`](./js/package.json).

Then, when you push or merge the code with the incremented versions to master, the GitHub Action `deploy.yml` workflow will build and push the packages to PyPI and NPM.


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
