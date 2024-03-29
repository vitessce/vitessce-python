import argparse
from anndata import read_h5ad
from vitessce.data_utils import (
    to_uint8,
    sort_var_axis,
    optimize_adata,
)


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)

    # Store an expression matrix with only the highly variable genes.
    adata = adata[:, adata.var['highly_variable']].copy()

    # Reorder the genes axis after hierarchical clustering.
    leaf_list = sort_var_axis(adata.X, adata.var.index.values)
    adata = adata[:, leaf_list].copy()

    # Store expression matrix as uint8.
    adata.layers['X_uint8'] = to_uint8(adata.X, norm_along="var")

    adata = optimize_adata(
        adata,
        obs_cols=["CellType"],
        obsm_keys=["X_umap"],
        layer_keys=["X_uint8"]
    )

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
