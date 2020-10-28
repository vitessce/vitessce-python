Getting Started
################

Python `Jupyter notebooks <https://jupyter.org>`_ are an excellent way to
experiment with data science and visualization. Using the vitessce-jupyter
extension, you can use Vitessce directly from within a Jupyter notebook.

Installation
-------------

To use higlass within a Jupyter notebook you need to install a few packages
and enable the Jupyter extension:


.. code-block:: bash

    pip install jupyter
    pip install vitessce

    jupyter nbextension enable --py --sys-prefix vitessce

If you use `JupyterLab <https://jupyterlab.readthedocs.io/en/stable/>`_ you also have to run

.. code-block:: bash

    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install vitessce


Simplest Use Case
------------------

The simplest way to instantiate a Vitessce widget is to create a view config and pass the view config as a parameter to the widget constructor:

.. code-block:: python

  from vitessce import VitessceConfig, VitessceWidget

  vc = VitessceConfig()
  vw = VitessceWidget(vc)
  vw
