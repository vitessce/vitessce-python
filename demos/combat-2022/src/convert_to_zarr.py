import argparse
from anndata import read_h5ad
import numpy as np
from scipy import sparse
from vitessce.data_utils import to_uint8


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)

    adata.var['is_gene_expression'] = adata.var['feature_types'] == 'Gene Expression'
    adata.var['is_antibody_capture'] = adata.var['feature_types'] == 'Antibody Capture'

    # TODO: automate conversion to csc in optimize_adata function
    adata.obsm['X_gene_expression_uint8'] = to_uint8(adata[:, adata.var['is_gene_expression']].X, norm_along="global")
    if isinstance(adata.obsm['X_gene_expression_uint8'], sparse.spmatrix):
        adata.obsm['X_gene_expression_uint8'] = adata.obsm['X_gene_expression_uint8'].tocsc()
    adata.obsm['X_antibody_capture_uint8'] = to_uint8(adata[:, adata.var['is_antibody_capture']].X, norm_along="global")
    if isinstance(adata.obsm['X_antibody_capture_uint8'], sparse.spmatrix):
        adata.obsm['X_antibody_capture_uint8'] = adata.obsm['X_antibody_capture_uint8'].tocsc()

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
