import argparse
from anndata import read_h5ad
import numpy as np


def to_uint8(arr):
    arr *= 255.0 / arr.max()
    arr = arr.astype(np.dtype('uint8')).todense()
    return arr


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)
    
    adata.X = to_uint8(adata.X)

    adata.write_zarr(output_path)


if __name__ == '__main__':
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
    convert_h5ad_to_zarr(
        args.input,
        args.output
    )