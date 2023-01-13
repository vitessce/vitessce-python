import argparse
from anndata import read_h5ad
import numpy as np
import scanpy as sc
from scipy import sparse
from vitessce.data_utils import (
    to_diamond,
    to_uint8,
    optimize_adata,
)


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var['feature_name'].str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_hvg, max_value=3)

    adata.obsm['X_hvg'] = adata_hvg.X
    adata.obsm['X_hvg_uint8'] = to_uint8(adata_hvg.X, norm_along="var")

    num_cells = adata.obs.shape[0]
    adata.obsm['X_spatial'] = adata.obsm['X_spatial']
    adata.obsm['X_segmentations'] = np.zeros((num_cells, 4, 2))
    radius = 10
    for i in range(num_cells):
        adata.obsm['X_segmentations'][i, :, :] = to_diamond(adata.obsm['X_spatial'][i, 0], adata.obsm['X_spatial'][i, 1], radius)

    adata = optimize_adata(
        adata,
        obs_cols=["cell_type"],
        var_cols=["feature_name"],
        obsm_keys=["X_hvg", "X_hvg_uint8", "X_umap", "X_spatial", "X_segmentations"],
        layer_keys=[],
    )

    adata.write_zarr(output_path, chunks=[adata.shape[0], 10])


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
