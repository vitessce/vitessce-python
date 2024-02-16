import argparse
from anndata import read_h5ad
import scipy
import numpy as np
import pandas as pd
import sparse
import dask
import dask.array as da
from dask.distributed import Client, LocalCluster, progress
import platform
import os

def convert_h5ad_to_zarr(input_path, output_path, client):
    adata = read_h5ad(input_path)

    # Clear X so that we can write it ourselves manually
    X = adata.X #.copy()
    
    adata.X = None
    adata.write_zarr(output_path)

    assert isinstance(X, scipy.sparse.spmatrix)

    # Use dask to write the X matrix as a dense matrix
    X_sparse = sparse.GCXS.from_scipy_sparse(X)
    print("X_sparse")
    X_dask = da.from_array(X_sparse, chunks=(3000, 1))
    print("X_dask")
    X_dask = X_dask.map_blocks(lambda block: block.todense())
    print("map_blocks")
    X_delayed = X_dask.to_zarr(url=output_path, component="/X", overwrite=True, compute=False)
    print("to_zarr")
    X_future = client.persist(X_delayed)

    progress(X_future)


if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Input H5AD file'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='Output Zarr store'
    )
    args = parser.parse_args()

    # Dask
    # Check if this is running on O2
    IS_O2 = (platform.system() == "Linux")

    if IS_O2:
        O2_USER = os.environ["USER"]
        DASK_TEMP_DIR = f"/n/scratch/users/{O2_USER[0]}/{O2_USER}/vitessce-python-temp"
        os.makedirs(DASK_TEMP_DIR, exist_ok=True)
        dask.config.set({ "temporary_directory": DASK_TEMP_DIR })

    # Should request at least 96GB of memory for this job.
    cluster = LocalCluster(n_workers=2, threads_per_worker=2, memory_limit='4GB')
    client = Client(cluster)

    convert_h5ad_to_zarr(
        args.input,
        args.output,
        client,
    )
