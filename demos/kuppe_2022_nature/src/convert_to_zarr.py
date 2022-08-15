import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import read_h5ad, AnnData
import imageio
import zarr
import ome_zarr


def to_uint8(arr):
    arr *= 255.0 / arr.max()
    arr = arr.astype(np.dtype('uint8')).todense()
    return arr


def process_h5ad_files(args):

    visium_img = imageio.imread(args.visium_img)
    visium_df = pd.read_csv(args.visium_csv)

    # TODO: write visium_img to OME-Zarr
    # https://github.com/vitessce/vitessceR/blob/main/R/data_to_zarr.R#L146
    # TODO: add the tissue_positions_list columns to visium_adata.obs
    # TODO: use scale factors from scalefactors_json.json in the OME-Zarr?

    rna_adata = read_h5ad(args.input_rna)
    atac_adata = read_h5ad(args.input_atac)
    visium_adata = read_h5ad(args.input_visium_adata)

    rna_adata.X = to_uint8(rna_adata.X)
    atac_adata.X = to_uint8(atac_adata.X)

    rna_adata.obsm['X_umap'] = rna_adata.obsm['X_umap'].astype('<f4')
    atac_adata.obsm['X_umap'] = atac_adata.obsm['X_umap'].astype('<f4')

    joint_cols = ['cell_type', 'development_stage', 'disease', 'sex']
    joint_obs_df = pd.concat([
        rna_adata.obs[joint_cols],
        atac_adata.obs[joint_cols],
        visium_adata.obs[joint_cols]
    ])
    
    joint_adata = AnnData(obs=joint_obs_df)
    joint_adata.write_zarr(args.output_joint)

    rna_adata.write_zarr(args.output_rna)
    atac_adata.write_zarr(args.output_atac)

    # Visium processing
    num_cells = visium_adata.obs.shape[0]
    visium_adata.obsm['X_spatial'] = visium_adata.obsm['X_spatial'].astype('<f4')

    # Create segmentations
    def to_diamond(x, y, r):
        return np.array([[x, y + r], [x + r, y], [x, y - r], [x - r, y]])
    visium_adata.obsm['segmentations'] = np.zeros((num_cells, 4, 2), dtype=np.dtype('<f4'))
    radius = 50
    for i in range(num_cells):
        visium_adata.obsm['segmentations'][i, :, :] = to_diamond(visium_adata.obsm['X_spatial'][i, 0], visium_adata.obsm['X_spatial'][i, 1], radius)

    visium_adata.write_zarr(args.output_visium)



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
        '-iva',
        '--input_visium_adata',
        type=str,
        required=True,
        help='Input visium h5ad file'
    )
    parser.add_argument(
        '-ivi',
        '--input_visium_img',
        type=str,
        required=True,
        help='Input visium PNG file'
    )
    parser.add_argument(
        '-ivc',
        '--input_visium_csv',
        type=str,
        required=True,
        help='Input visium CSV file'
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
    parser.add_argument(
        '-ov',
        '--output_visium',
        type=str,
        required=True,
        help='Output visium h5ad file'
    )
    args = parser.parse_args()
    process_h5ad_files(args)
