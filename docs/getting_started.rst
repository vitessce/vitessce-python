Getting Started
################

The ``vitessce`` Python package has two main parts:

* `Widget API <widget_api.html>`_: use Vitessce directly from with Jupyter Notebook or Jupyter Lab
* `Config API <config_api.html>`_: create and edit Vitessce view configs using Python syntax

.. image:: jupyter_screenshot.png

Installation
-------------

To use Vitessce within a Jupyter notebook you need to install a few packages
and enable the Jupyter extension:


.. code-block:: bash
    :emphasize-lines: 1

    pip install vitessce
    pip install jupyter
    jupyter nbextension enable --py --sys-prefix vitessce

If you use `JupyterLab <https://jupyterlab.readthedocs.io/en/stable/>`_ you also have to run

.. code-block:: bash

    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install vitessce


Simplest Widget Use Case
------------------

The simplest way to instantiate a Vitessce widget is to create a view config based on a single-cell dataset object (from which data types and visualization types can be inferred automatically) and pass the view config instance as a parameter to the widget constructor:

.. code-block:: python

  from vitessce import VitessceConfig, VitessceWidget

  vc = VitessceConfig.from_object(my_scanpy_object)
  vw = VitessceWidget(vc)
  vw

To customize the view config passed into the widget (or to define a view config manually), please see the documentation for the ``VitessceConfig`` API.