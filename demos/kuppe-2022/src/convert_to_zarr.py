import argparse
import numpy as np
import pandas as pd
import json
from anndata import read_h5ad, AnnData
from scipy import sparse
import imageio.v2 as imageio
from vitessce.data_utils import (
    to_diamond,
    to_uint8,
    rgb_img_to_ome_zarr,
    optimize_adata,
)


def process_h5ad_files(args):
    img_arr = imageio.imread(args.input_visium_img)
    img_arr = np.transpose(img_arr, axes=(2, 1, 0))  # xyc to cyx

    visium_df = pd.read_csv(args.input_visium_csv, header=None)
    visium_df = visium_df.rename(columns={
        0: "spot_id",
        4: "X",
        5: "Y"
    })

    visium_df = visium_df.set_index('spot_id')

    # Write img_arr to OME-Zarr
    # https://github.com/vitessce/vitessceR/blob/main/R/data_to_zarr.R#L146
    rgb_img_to_ome_zarr(
        img_arr,
        args.output_visium_ome,
        img_name="GT_IZ_P9",
        axes="cyx",
        chunks=(1, 256, 256)
    )

    visium_adata = read_h5ad(args.input_visium_adata)

    # Add the tissue_positions_list columns to visium_adata.obs
    visium_adata.obs['X'] = visium_adata.obs.apply(lambda row: visium_df.at[row.name, 'X'], axis='columns')
    visium_adata.obs['Y'] = visium_adata.obs.apply(lambda row: visium_df.at[row.name, 'Y'], axis='columns')

    rna_adata = read_h5ad(args.input_rna)
    atac_adata = read_h5ad(args.input_atac)

    rna_adata.layers['X_uint8'] = to_uint8(rna_adata.X, norm_along="global")
    atac_adata.layers['X_uint8'] = to_uint8(atac_adata.X, norm_along="global")

    visium_adata.layers['X_uint8'] = to_uint8(visium_adata.X, norm_along="var")

    joint_cols = ['cell_type', 'development_stage', 'disease', 'sex']
    joint_obs_df = pd.concat([
        rna_adata.obs[joint_cols],
        atac_adata.obs[joint_cols],
        visium_adata.obs[joint_cols]
    ])

    joint_adata = AnnData(obs=joint_obs_df)
    joint_adata.write_zarr(args.output_joint)

    joint_adata = optimize_adata(
        joint_adata,
        obs_cols=["cell_type", "development_stage", "disease", "sex"]
    )
    rna_adata = optimize_adata(
        rna_adata,
        obsm_keys=["X_umap", "X_pca"],
        var_cols=["feature_name"],
        layer_keys=["X_uint8"],
        optimize_X=False,
    )

    atac_adata = optimize_adata(
        atac_adata,
        obsm_keys=["X_umap"],
        var_cols=["feature_name"],
        layer_keys=["X_uint8"],
        optimize_X=False,
    )

    # Use chunks in case data is not sparse.
    rna_adata.write_zarr(args.output_rna, [rna_adata.shape[0], 10])
    atac_adata.write_zarr(args.output_atac, [atac_adata.shape[0], 10])

    # Visium processing
    num_cells = visium_adata.obs.shape[0]
    visium_adata.obsm['spatial'] = visium_adata.obsm['X_spatial']
    visium_adata.obsm['xy'] = visium_adata.obs[['X', 'Y']].values

    with open(args.input_visium_img_scalefactors) as f:
        scale_factor = json.load(f)['tissue_hires_scalef']
    visium_adata.obsm['xy_scaled'] = visium_adata.obsm['xy'] * scale_factor

    # Create segmentations
    visium_adata.obsm['segmentations'] = np.zeros((num_cells, 4, 2))
    visium_adata.obsm['xy_segmentations'] = np.zeros((num_cells, 4, 2))
    visium_adata.obsm['xy_segmentations_scaled'] = np.zeros((num_cells, 4, 2))
    radius = 35
    for i in range(num_cells):
        visium_adata.obsm['segmentations'][i, :, :] = to_diamond(visium_adata.obsm['spatial'][i, 0], visium_adata.obsm['spatial'][i, 1], radius)
        visium_adata.obsm['xy_segmentations'][i, :, :] = to_diamond(visium_adata.obsm['xy'][i, 0], visium_adata.obsm['xy'][i, 1], radius)
        visium_adata.obsm['xy_segmentations_scaled'][i, :, :] = to_diamond(visium_adata.obsm['xy_scaled'][i, 0], visium_adata.obsm['xy_scaled'][i, 1], radius * scale_factor)

    visium_adata = optimize_adata(
        visium_adata,
        obs_cols=["cell_type", "development_stage", "disease", "sex"],
        obsm_keys=["xy_scaled", "xy_segmentations_scaled"],
        var_cols=["feature_name"],
        layer_keys=["X_uint8"],
    )
    visium_adata.write_zarr(args.output_visium_adata)


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
        '-ivs',
        '--input_visium_img_scalefactors',
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
        '-ova',
        '--output_visium_adata',
        type=str,
        required=True,
        help='Output visium h5ad zarr store'
    )
    parser.add_argument(
        '-ovo',
        '--output_visium_ome',
        type=str,
        required=True,
        help='Output visium OME zarr store'
    )
    args = parser.parse_args()
    process_h5ad_files(args)
