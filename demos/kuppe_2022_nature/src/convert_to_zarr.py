import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import read_h5ad, AnnData


def to_uint8(arr):
    arr *= 255.0 / arr.max()
    arr = arr.astype(np.dtype('uint8')).todense()
    return arr


def process_h5ad_files(args):
    rna_adata = read_h5ad(args.input_rna)
    atac_adata = read_h5ad(args.input_atac)

    rna_adata.X = to_uint8(rna_adata.X)
    atac_adata.X = to_uint8(atac_adata.X)

    joint_cols = ['cell_type', 'development_stage', 'disease', 'sex']
    joint_obs_df = pd.concat([ rna_adata.obs[joint_cols], atac_adata.obs[joint_cols] ])
    
    joint_adata = AnnData(obs=joint_obs_df)
    joint_adata.write_zarr(args.output_joint)

    rna_adata.write_zarr(args.output_rna)
    atac_adata.write_zarr(args.output_atac)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-ir',
        '--input_rna',
        type=str,
        required=True,
        help='Input snRNA-seq h5ad file'
    )
    parser.add_argument(
        '-ia',
        '--input_atac',
        type=str,
        required=True,
        help='Input snATAC-seq h5ad file'
    )
    parser.add_argument(
        '-or',
        '--output_rna',
        type=str,
        required=True,
        help='Output snRNA-seq zarr store'
    )
    parser.add_argument(
        '-oa',
        '--output_atac',
        type=str,
        required=True,
        help='Output snATAC-seq zarr store'
    )
    parser.add_argument(
        '-oj',
        '--output_joint',
        type=str,
        required=True,
        help='Output joint RNA+ATAC zarr store'
    )
    args = parser.parse_args()
    process_h5ad_files(args)
