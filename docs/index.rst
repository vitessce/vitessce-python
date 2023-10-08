Vitessce is a visual integration tool for exploration of spatial single-cell experiments.
To learn more about the features of Vitessce, please visit our `core docs <http://vitessce.io>`_.

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

Installation requires Python 3.9 or greater. Install from `PyPI <https://pypi.org/project/vitessce>`_:

.. code-block:: bash

    pip install vitessce[all]


Notes:

* The most minimal ``pip install vitessce`` installs only those dependencies required to use the configuration and wrapper classes, such as ``VitessceConfig`` and ``AnnDataWrapper``.
* The second-most minimal ``pip install vitessce[all]`` additionally installs the dependencies required to use the Jupyter widget, enabling ``VitessceConfig.widget()``.
* Data conversion dependencies for usage of ``vitessce.data_utils`` submodules must be installed explicitly as described below.

To use functions from ``vitessce.data_utils.{submodule}``, the name of each submodule intended to be used must be specified in `extras <https://peps.python.org/pep-0508/#extras>`_:

.. code-block:: bash

    pip install vitessce[anndata,multivec,ome_tiff,ome_zarr,ucsc_cellbrowser]
  
For example, to install the Jupyter widget and the data conversion dependencies for ``vitessce.data_utils.anndata``:

.. code-block:: bash

    pip install vitessce[all,anndata]




Widget Compatibility
--------------------

The Vitessce widget is compatible with the following interactive Python platforms:

* JupyterLab ``>=3.0.0``
* Jupyter Notebook (classic) ``>=1.0.0``


.. toctree::
   :maxdepth: 2
   :hidden:

   self
   widget_examples
   data_examples
   api_config
   api_data
   api_data_utils
   data_options
   screenshots



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
