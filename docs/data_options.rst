Data location options
#####################

There are several possibilities for how to configure data for the Vitessce Jupyter widget depending on:

* where the Python process powering Jupyter is running (i.e., on which machine, relative to the web browser in which Jupyter is being accessed)
* where the data is relative to the Python process for Jupyter
* where the web browser accessing Jupyter is relative to the Python process for Jupyter

The following sections describe the various possibilities and the implications that each has on Vitessce's automatic data export functionality.
The goal of automatic data export is: given an input Vitessce configuration (i.e., ``VitessceConfig`` instance), to generate a directory of data files that, when served via HTTP, can be accessed by Vitessce in a web browser using the corresonding JSON-serialized configuration (i.e., ``VitessceConfig.to_dict()``).

====================================
Jupyter process: local; Files: local
====================================

In this case, you ran ``jupyter lab`` in a terminal local to your computer (i.e., not a cluster or remote machine) and the data files you want to visualize are located on the same machine.

-------------------------------------------------------
Configuration via file paths relative to a ``base_dir``
-------------------------------------------------------

For more information, see the VitessceConfig `constructor <api_config.html#vitessce.config.VitessceConfig>`_ or the `Configure relative to a base_dir <notebooks/widget_brain_with_base_dir.html>`_ example notebook.

**Note**: Export (and thereby copying of files) is not required, as ``base_dir`` is already equivalent to an exported data directory, which the file URLs in the configuration will be relative to.

-------------------------------------
Configuration via relative file paths
-------------------------------------

**Note**: Export requires copying, as we need to generate a single data directory for all of the files, which the file URLs in the configuration will be relative to.

-------------------------------------
Configuration via absolute file paths
-------------------------------------

**Note**: Export requires copying, as we need to generate a single data directory for all of the files, which the file URLs in the configuration will be relative to.


==========================================================
Jupyter process: anywhere; Files: remote & served via HTTP
==========================================================

In this case, Jupyter may be running anywhere (e.g., local, or remote machine, or service like Colab/Binder).
Files are already being served via HTTP on internet-accessible URLs (i.e., non-localhost).

For more information about how to host data for Vitessce on remote servers, please see http://vitessce.io/docs/data-hosting/.


**Note**: Exporting the data is not required (nor possible) as the files are already referenced via HTTP URLs in the configuration.

--------------
Range requests
--------------

Certain file formats (e.g., OME-TIFF) are loaded by Vitessce via HTTP range requests, so ensure your hosting service supports these when using such file formats.

========================================================
Jupyter process: local; Files: remote & accessed via SSH
========================================================

In this case, you can serve the files using a local HTTP server such as `http-server <https://github.com/http-party/http-server>`_ on a particular port, and then use port forwarding when initializing the SSH session.

If running the HTTP server within an interactive compute session (e.g., with SLURM), you may also need to configure port forwarding when initializing the compute session.

.. code-block:: bash

    # Replace REMOTE_PORT and LOCAL_PORT below with the desired port numbers.
    # Replace my_username and cluster.university.edu.
    ssh -L REMOTE_PORT:127.0.0.1:LOCAL_PORT my_username@cluster.university.edu
    # Optionally SSH to a particular login node first.
    ssh -L REMOTE_PORT:127.0.0.1:REMOTE_PORT login01
    # Optionally start an interactive compute session first, but make sure tunneling/port forwarding is enabled.
    # Note: this command may be different depending on the cluster and job management system.
    srun -t 0-3:00 --pty -p interactive --tunnel REMOTE_PORT:REMOTE_PORT /bin/bash
    # cd to some directory with files to serve.
    http-server --cors='*' --port REMOTE_PORT .


Then, you can configure Vitessce using localhost HTTP urls as if the files are being served locally (i.e., on ``http://localhost:LOCAL_PORT/``).

**Note**: Exporting the data is not possible as the files are already referenced via HTTP URLs in the configuration.

=====================================================================================
Jupyter process: remote & accessed via SSH; Files: on same machine as Jupyter process
=====================================================================================

This case can be treated almost the same as when the Jupyter process is local and the files are local (the first case).

However, when accessing the notebook from your local web browser, you may need to use the ``proxy`` parameter:

.. code-block:: python

    vc.widget(proxy=True)


===============================================================
Jupyter process: remote service like Colab/Binder; Files: local
===============================================================

Unfortunately, this will not work because the remote server cannot access your local files.

===================================================================================
Jupyter process: remote service like Colab/Binder; Files: remote & accessed via SSH
===================================================================================

Unfortunately, this will not work because the remote server cannot access the files that are on another machine behind SSH.

========================================================================
Jupyter process: anywhere; Files: anywhere that can be accessed via Zarr
========================================================================

If the data is readable via Zarr (i.e., `zarr.storage.*Store`) and the Jupyter process can access the store contents, then the Vitessce widget can access the data by specifying the Zarr store as the data source for Vitessce data wrapper class instances.
This is currently supported for the ``AnnDataWrapper`` class using the ``adata_store`` parameter (as opposed to ``adata_path`` or ``adata_url``).

.. code-block:: python

    from vitessce import VitessceConfig, AnnDataWrapper

    # ...
    adata.write_zarr("my_store.adata.zarr")

    vc = VitessceConfig(name="My Vitessce Configuration")
    vc.add_dataset(name="My Dataset").add_object(AnnDataWrapper(
        adata_store="my_store.adata.zarr",
        # ...
    ))
    # ...
    vc.widget()


Or, with a Zarr store instance (instead of a local string path to a DirectoryStore):

.. code-block:: python

    import zarr
    from vitessce import VitessceConfig, AnnDataWrapper

    # ...
    store = zarr.storage.FSStore("s3://my_bucket/path/to/my_store.adata.zarr")

    vc = VitessceConfig(name="My Vitessce Configuration")
    vc.add_dataset(name="My Dataset").add_object(AnnDataWrapper(
        adata_store=store,
        # ...
    ))
    # ...
    vc.widget()

