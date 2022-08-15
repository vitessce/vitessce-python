import argparse
import scanpy as sc
import numpy as np
from anndata import read_h5ad


def process_h5ad_files(args):
    rna_adata = read_h5ad(args.input_rna)
    atac_adata = read_h5ad(args.input_atac)

    print(rna_adata)
    print(atac_adata)

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
    args = parser.parse_args()
    process_h5ad_files(args)
