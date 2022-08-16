import argparse
import numpy as np
import pandas as pd
from anndata import read_h5ad, AnnData
import imageio.v2 as imageio
import zarr
from ome_zarr.writer import write_image


def to_uint8(arr):
    arr *= 255.0 / arr.max()
    arr = arr.astype(np.dtype('uint8')).todense()
    return arr


def to_uint8_norm(arr):
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

    default_window = {
        "start": 0,
        "min": 0,
        "max": 255,
        "end": 255
    }

    z_root = zarr.open_group(args.output_visium_ome)
    write_image(
        image=img_arr,
        group=z_root,
        axes="cyx",
        omero={
            "name": "GT_IZ_P9",
            "version": "0.3",
            "rdefs": {},
            "channels": [
                {
                    "label": "R",
                    "color": "FF0000",
                    "window": default_window
                },
                {
                    "label": "G",
                    "color": "00FF00",
                    "window": default_window
                },
                {
                    "label": "B",
                    "color": "0000FF",
                    "window": default_window
                }
            ]
        },
        chunks=(1, 256, 256)
    )

    visium_adata = read_h5ad(args.input_visium_adata)

    # Add the tissue_positions_list columns to visium_adata.obs
    visium_adata.obs['X'] = visium_adata.obs.apply(lambda row: visium_df.at[row.name, 'X'], axis='columns')
    visium_adata.obs['Y'] = visium_adata.obs.apply(lambda row: visium_df.at[row.name, 'Y'], axis='columns')

    # TODO: use scale factors from scalefactors_json.json in the OME-Zarr?

    rna_adata = read_h5ad(args.input_rna)
    atac_adata = read_h5ad(args.input_atac)

    rna_adata.X = to_uint8(rna_adata.X)
    atac_adata.X = to_uint8(atac_adata.X)
    visium_adata.X = to_uint8_norm(visium_adata.X.todense())

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
    visium_adata.obsm['spatial'] = visium_adata.obsm['X_spatial'].astype('<f4')
    visium_adata.obsm['xy'] = visium_adata.obs[['X', 'Y']].values.astype('<f4')

    scale_factor = 0.20319009
    visium_adata.obsm['xy_scaled'] = visium_adata.obsm['xy'] * scale_factor

    # Create segmentations
    def to_diamond(x, y, r):
        return np.array([[x, y + r], [x + r, y], [x, y - r], [x - r, y]])
    visium_adata.obsm['segmentations'] = np.zeros((num_cells, 4, 2), dtype=np.dtype('<f4'))
    visium_adata.obsm['xy_segmentations'] = np.zeros((num_cells, 4, 2), dtype=np.dtype('<f4'))
    visium_adata.obsm['xy_segmentations_scaled'] = np.zeros((num_cells, 4, 2), dtype=np.dtype('<f4'))
    radius = 35
    for i in range(num_cells):
        visium_adata.obsm['segmentations'][i, :, :] = to_diamond(visium_adata.obsm['spatial'][i, 0], visium_adata.obsm['spatial'][i, 1], radius)
        visium_adata.obsm['xy_segmentations'][i, :, :] = to_diamond(visium_adata.obsm['xy'][i, 0], visium_adata.obsm['xy'][i, 1], radius)
        visium_adata.obsm['xy_segmentations_scaled'][i, :, :] = to_diamond(visium_adata.obsm['xy_scaled'][i, 0], visium_adata.obsm['xy_scaled'][i, 1], radius * scale_factor)

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
