vitessce-jupyter
===============================

Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce

Installation
------------

To install use pip:

    $ pip install vitessce
    $ jupyter nbextension enable --py --sys-prefix vitessce

To install for jupyterlab

    $ jupyter labextension install vitessce

For a development installation (requires npm),

    $ git clone https://github.com/keller-mark/vitessce-jupyter.git
    $ cd vitessce-jupyter
    $ pip install -e .
    $ jupyter nbextension install --py --symlink --sys-prefix vitessce
    $ jupyter nbextension enable --py --sys-prefix vitessce
    $ cd js
    $ jupyter labextension install

When actively developing your extension, build Jupyter Lab with the command:

    $ jupyter lab --watch

This takes a minute or so to get started, but then automatically rebuilds JupyterLab when your javascript changes.

Note on first `jupyter lab --watch`, you may need to touch a file to get Jupyter Lab to open.

