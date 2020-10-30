# vitessce-jupyter

Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce

🚧 work in progress 👷


## Installation

To install use pip:

    $ pip install vitessce
    $ jupyter nbextension enable --py --sys-prefix vitessce

To install for jupyterlab

    $ jupyter labextension install @jupyter-widgets/jupyterlab-manager
    $ jupyter labextension install vitessce


## Getting started

Explore our [example notebooks](./notebooks/).
These contain demos of different use cases and integrations with single-cell data analysis packages.


## Development

For a development installation (requires npm),

    $ git clone https://github.com/keller-mark/vitessce-jupyter.git
    $ cd vitessce-jupyter
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
python -m unittest
```


## Documentation

```sh
make html
```


## Resources

- [ipywidget docs: Building a Custom Widget](https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20Custom.html)
- [ipywidget docs: Low Level Widget Tutorial](https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Low%20Level.html)
- [ipywidget example: ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet)
- [ipywidget example (with React): ipymaterialui](https://github.com/maartenbreddels/ipymaterialui)
- [ipywidget example (with React): higlass-python](https://github.com/higlass/higlass-python)
- [ipywidget cookiecutter](https://github.com/jupyter-widgets/widget-cookiecutter)
- [Sphinx: Getting Started](https://www.sphinx-doc.org/en/master/usage/quickstart.html)
- [Read the Docs Sphinx Theme](https://github.com/readthedocs/sphinx_rtd_theme)