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

Installation requires Python 3.8 or greater. Install from `PyPI <https://pypi.org/project/vitessce>`_:

.. code-block:: bash

    pip install vitessce


Widget Compatibility
--------------------

The Vitessce widget is compatible with the following interactive Python platforms:

* JupyterLab ``>=3.0.0``
* Jupyter Notebook (classic) ``>=1.0.0``


Installation with optional dependencies
---------------------

To use the widget through a proxy (e.g. a cloud notebook service like Binder). See widget parameter ``proxy``.

.. code-block:: bash

    pip install vitessce[proxy]


.. toctree::
   :maxdepth: 2
   :hidden:

   self
   widget_examples
   data_examples
   api_config
   api_data
   screenshots



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
