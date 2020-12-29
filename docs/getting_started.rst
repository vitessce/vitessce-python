Getting Started
################

The ``vitessce`` Python package includes:

* **Vitessce as a Jupyter Widget**

  * Use Vitessce directly within Jupyter Notebook or Jupyter Lab as an interactive widget

* **View config API**

  * Create and edit Vitessce configurations using Python object-oriented syntax

* **Data preparation**

  * Use our data conversion wrapper classes to process data stored in common single-cell file types including AnnData and SnapATAC.


Installation
-------------

To use Vitessce within a Jupyter notebook you need to install a few packages
and enable the Jupyter extension:


.. code-block:: bash

    pip install vitessce
    jupyter nbextension enable --py --sys-prefix vitessce

To use the widget in `JupyterLab <https://jupyterlab.readthedocs.io/en/stable/>`_ also run the following:

.. code-block:: bash

    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install vitessce-jupyter


Optional dependencies
---------------------

The following dependencies are optional, and only required when using certain data wrapper classes.

* ``anndata>=0.7.4`` - Required for AnnData support with the :class:`~vitessce.wrappers.AnnDataWrapper` class.
* ``loompy>=3.0.6`` - Required for Loom support with the :class:`~vitessce.wrappers.LoomWrapper` class.

