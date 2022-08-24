import argparse
from anndata import read_h5ad
import numpy as np
import pandas as pd
import scipy.cluster


def clean_gexp(adata):
    adata.X = adata.X.todense()
    gexp_arr = adata.X
    gexp_df = adata.to_df()

    # Re-scale the gene expression values between 0 and 255 (one byte ints).
    gexp_arr_min = gexp_arr.min()
    gexp_arr_max = gexp_arr.max()
    gexp_arr_range = gexp_arr_max - gexp_arr_min
    gexp_arr_ratio = 255 / gexp_arr_range
    gexp_norm_arr = (gexp_arr - gexp_arr_min) * gexp_arr_ratio

    # Perform hierarchical clustering along the genes axis.
    Z = scipy.cluster.hierarchy.linkage(gexp_norm_arr.T, method="ward")
    labels = adata.var.index.values

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = labels[leaf_index_list].tolist()

    # Create a new *ordered* gene expression dataframe.
    gexp_norm_df = pd.DataFrame(
        index=gexp_df.index.values.tolist(),
        columns=gexp_df.columns.values.tolist(),
        data=gexp_norm_arr
    )
    adata.X = gexp_norm_df.values

    # Sort by selecting along the gene axis of the AnnData object
    adata_sorted = adata[:, leaf_list].copy()
    adata.X = adata.X.astype(np.dtype('uint8'))
    return adata_sorted


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)

    # Store an expression matrix with only the highly variable genes.
    adata = adata[:, adata.var['highly_variable']].copy()

    adata = clean_gexp(adata)

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
