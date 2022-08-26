import argparse
from anndata import read_h5ad
import numpy as np
from vitessce.data_utils import to_uint8, to_dense, optimize_adata


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path, backed="r+")

    adata.var['is_gene_expression'] = adata.var['feature_types'] == 'Gene Expression'
    adata.var['is_antibody_capture'] = adata.var['feature_types'] == 'Antibody Capture'

    adata.obsm['X_gene_expression'] = adata.X[:, adata.var['is_gene_expression']]
    adata.obsm['X_antibody_capture'] = adata.X[:, adata.var['is_antibody_capture']]
    adata.obsm['X_gene_expression_uint8'] = to_uint8(adata.obsm['X_gene_expression'], norm_along="global")
    adata.obsm['X_antibody_capture_uint8'] = to_uint8(adata.obsm['X_antibody_capture'], norm_along="global")

    print("running optimize_adata")
    adata = optimize_adata(
        adata,
        ignore_X=True,
        obs_cols=[
            "cell_type", "cluster",
            "major_subset", "minor_subset",
            "COMBAT_ID", "scRNASeq_sample_ID",
            "GEX_region",
            "Source", "Sex", "Age", "Hospitalstay", "TimeSinceOnset"
        ],
        obsm_keys=[
            "X_gene_expression", "X_antibody_capture",
            "X_gene_expression_uint8", "X_antibody_capture_uint8",
            "X_umap",
        ],
        var_cols=["feature_types", "is_gene_expression", "is_antibody_capture", "gene_ids"],
        layer_keys=[],
    )
    print("done optimize_adata")

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
