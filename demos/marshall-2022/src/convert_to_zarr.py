import argparse
from anndata import read_h5ad
import numpy as np
import scanpy as sc


def to_uint8(arr):
    num_cells = arr.shape[0]
    min_along_genes = arr.min(axis=0)
    max_along_genes = arr.max(axis=0)
    range_per_gene = max_along_genes - min_along_genes
    ratio_per_gene = 255.0 / range_per_gene

    norm_along_genes_arr = np.multiply(
        (arr - np.tile(min_along_genes, (num_cells, 1))),
        np.tile(ratio_per_gene, (num_cells, 1))
    )
    return norm_along_genes_arr.astype(np.dtype('uint8'))


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)
    print(adata)

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

    adata.obsm['X_hvg'] = to_uint8(adata_hvg.X)

    def to_diamond(x, y, r):
        return np.array([[x, y + r], [x + r, y], [x, y - r], [x - r, y]])

    num_cells = adata.obs.shape[0]
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'].astype(np.dtype('uint16'))
    adata.obsm['X_segmentations'] = np.zeros((num_cells, 4, 2), dtype=np.dtype('uint16'))
    radius = 10
    for i in range(num_cells):
        adata.obsm['X_segmentations'][i, :, :] = to_diamond(adata.obsm['X_spatial'][i, 0], adata.obsm['X_spatial'][i, 1], radius)

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
