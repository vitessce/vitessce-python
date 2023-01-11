import argparse
from anndata import read_zarr
import json
import pandas as pd


def convert_to_csv(args):
    adata = read_zarr(args.input)

    cells_df = adata.obs
    cells_df['TSNE_1'] = adata.obsm['X_tsne'][:, 0]
    cells_df['TSNE_2'] = adata.obsm['X_tsne'][:, 1]
    cells_df['PCA_1'] = adata.obsm['X_pca'][:, 0]
    cells_df['PCA_2'] = adata.obsm['X_pca'][:, 1]
    cells_df['UMAP_1'] = adata.obsm['X_umap'][:, 0]
    cells_df['UMAP_2'] = adata.obsm['X_umap'][:, 1]
    cells_df['X'] = adata.obsm['X_centroid'][:, 0]
    cells_df['Y'] = adata.obsm['X_centroid'][:, 1]
    cells_df.index = cells_df.index.rename("cell_id")

    segmentations = {}
    for i, cell_id in enumerate(adata.obs.index.values.tolist()):
        segmentations[cell_id] = [
            [int(coord) for coord in xy]
            for xy in adata.obsm["X_segmentations"][i, :, :]
        ]

    matrix_df = pd.DataFrame(
        index=adata.obs.index.values.tolist(),
        columns=adata.var.index.values.tolist(),
        data=adata.X
    )
    matrix_df.index = matrix_df.index.rename("cell_id")

    cells_df.to_csv(args.output_cells, index=True)
    matrix_df.to_csv(args.output_matrix, index=True)

    with open(args.output_segmentations, "w") as f:
        json.dump(segmentations, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input AnnData-Zarr store",
    )
    parser.add_argument(
        "-oc",
        "--output_cells",
        type=str,
        required=True,
        help="Output cells.csv file"
    )
    parser.add_argument(
        "-os",
        "--output_segmentations",
        type=str,
        required=True,
        help="Output obsSegmentations.json file",
    )
    parser.add_argument(
        "-om",
        "--output_matrix",
        type=str,
        required=True,
        help="Output obsFeatureMatrix.csv file",
    )
    args = parser.parse_args()
    convert_to_csv(args)
