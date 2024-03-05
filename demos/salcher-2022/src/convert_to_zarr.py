import argparse
from anndata import read_h5ad
import scipy
import numpy as np
import pandas as pd
import platform
import os
import zarr
import math
from os.path import join

def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)

    for adata_size in [10, 100, 1000, 10000, 100000, 1000000]:
        adata_small = adata[:adata_size, :].copy()
        adata_small.write_zarr(join(output_path, f"small_{adata_size}.h5ad.zarr"))
        adata_small.write_csvs(join(output_path, f"small_{adata_size}.csvs_dir"), skip_data=False)

        print(f"done for {adata_size}")

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
        help='Output directory'
    )
    args = parser.parse_args()

    convert_h5ad_to_zarr(
        args.input,
        args.output,
    )
