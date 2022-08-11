import argparse
import json
from anndata import read_h5ad
import numpy as np


def to_uint8(arr):
    arr *= 255.0/arr.max()
    arr = arr.astype(np.dtype('uint8')).todense()
    return arr

def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)
    
    adata.var['is_gene_expression'] = adata.var['feature_types'] == 'Gene Expression'
    adata.var['is_antibody_capture'] = adata.var['feature_types'] == 'Antibody Capture'

    adata.obsm['X_gene_expression'] = to_uint8(adata.X[:, adata.var['is_gene_expression']])
    adata.obsm['X_antibody_capture'] = to_uint8(adata.X[:, adata.var['is_antibody_capture']])

    adata.obsm['X_umap'] = adata.obsm['X_umap'].astype(np.dtype('<f4'))

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
