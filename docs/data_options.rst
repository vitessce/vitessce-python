Data location options
###############

There are several possibilities for how to configure data for the Vitessce Jupyter widget depending on:

* where the Python process powering Jupyter is running
* where the data is located relative to the Python process for Jupyter
* where the web browser accessing Jupyter is located relative to the Python process for Jupyter

=================
Jupyter process: local; Files: local
=================

In this case, you ran `jupyter lab` in a terminal local to your computer (i.e., not a cluster or remote machine) and the data files you want to visualize are located on the same machine.

-----------------
Configuration via relative file paths
-----------------

Export requires copying

-----------------
Configuration via absolute file paths
-----------------

Export requires copying

-----------------
Configuration via file paths relative to a `base_dir`
-----------------

Export not required; `base_dir` is equivalent to an exported data directory


=================
Jupyter process: anywhere; Files: remote & already being served via HTTP
=================

In this case, Jupyter is running anywhere (e.g., local, or remote machine, or service like Colab/Binder).
Files are being served via HTTP on internet-accessible URLs (i.e., non-localhost).

For more information about how to host data for Vitessce on remote servers, please see http://vitessce.io/docs/data-hosting/.


Export not required; http urls already absolute


-----------------
Range requests
-----------------

Certain file formats (e.g., OME-TIFF) are loaded by Vitessce via HTTP range requests, so ensure your hosting service supports these when using such file formats.

=================
Jupyter process: local; Files: remote & accessed via SSH
=================

In this case, you can serve the files using a local HTTP server such as `http-server` on a particular port, and then use port forwarding when initializing the SSH session.

If running the HTTP server within an interactive compute session (e.g., with SLURM), you may also need to configure port forwarding when initializing the compute session.

TODO: copy code snippet from https://github.com/keller-mark/snippets#serve-directory-of-files
TODO: link to docs to install `http-server`

Then, you can configure Vitessce using localhost HTTP urls as if the files are being served locally.

Note: Automatic export not possible.

=================
Jupyter process: remote & accessed via SSH; Files: on same machine as Jupyter process
=================

This case can be treated almost the same as when the Jupyter process is local and the files are local (the first case).

However, when accessing the notebook from your local web browser, you may need to use the `proxy` parameter:

```python
vc.widget(proxy=True)
```

  

=================
Jupyter process: remote service like Colab/Binder; Files: local
=================

Unfortunately, this will not work.

=================
Jupyter process: remote service like Colab/Binder; Files: remote & accessed via SSH
=================

Unfortunately, this will not work.



