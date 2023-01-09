import argparse
import json
from anndata import read_h5ad
import pandas as pd
from vitessce.data_utils import (
    to_uint8,
    sort_var_axis,
    optimize_adata,
)

from constants import COLUMNS
from generate_cell_sets import (
    generate_leiden_cluster_cell_sets,
    generate_cell_type_cell_sets
)
from utils import (
    merge_cell_sets_trees
)


def convert_to_zarr_and_json(
    input_cells_h5ad, input_annotations_csv, input_cl_obo_file,
    output_cells, output_cell_sets,
):

    adata = read_h5ad(input_cells_h5ad)
    adata.obs['leiden'] = adata.obs['leiden'].apply(lambda i: f"Cluster {str(i).zfill(2)}")

    annotation_df = pd.read_csv(input_annotations_csv)
    annotation_df = annotation_df.set_index(COLUMNS.CELL_ID.value)

    adata.obs[COLUMNS.ANNOTATION.value] = adata.obs.apply(lambda row: annotation_df.at[row.name, COLUMNS.ANNOTATION.value], axis='columns')
    adata.obs[COLUMNS.PREDICTION_SCORE.value] = adata.obs.apply(lambda row: annotation_df.at[row.name, COLUMNS.PREDICTION_SCORE.value], axis='columns')

    # Remove annotations with NaN prediction scores
    not_nan_cells = pd.notna(adata.obs[COLUMNS.PREDICTION_SCORE.value])
    adata = adata[not_nan_cells, :].copy()

    # Generate data for .cell_sets.json
    df = adata.obs
    df.index = df.index.rename(COLUMNS.CELL_ID.value)
    df = df.reset_index()

    leiden_cell_sets = generate_leiden_cluster_cell_sets(df)
    # The cell type annotations for this dataset are hierarchical.
    cell_type_cell_sets = generate_cell_type_cell_sets(
        df,
        input_cl_obo_file
    )

    # Merge the Leiden Cluster and Cell Type Annotation cell sets.
    cell_sets = merge_cell_sets_trees(
        leiden_cell_sets,
        cell_type_cell_sets
    )

    # These hierarchical cell type annotations must be represented
    # in a JSON format (rather than using a column-based representation
    # like AnnData.obs) because the tree does not have a uniform height.
    with open(output_cell_sets, 'w') as f:
        json.dump(cell_sets, f, indent=1)

    # Round prediction score for labels.
    adata.obs[COLUMNS.PREDICTION_SCORE.value] = adata.obs[COLUMNS.PREDICTION_SCORE.value].apply(lambda val: '{:.2g}'.format(val))

    # Reorder the genes axis after hierarchical clustering.
    leaf_list = sort_var_axis(adata.X, adata.var.index.values)
    adata = adata[:, leaf_list].copy()

    # Store expression matrix as uint8
    adata.layers['X_uint8'] = to_uint8(adata.X, norm_along="var")

    adata = optimize_adata(
        adata,
        obs_cols=["prediction_score"],
        obsm_keys=["X_umap"],
        layer_keys=["X_uint8"],
    )

    adata.write_zarr(output_cells)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-ic',
        '--input_cells_h5ad',
        type=str,
        required=True,
        help='Input h5ad file'
    )
    parser.add_argument(
        '-ia',
        '--input_annotations_csv',
        type=str,
        required=True,
        help='Input Arrow file'
    )
    parser.add_argument(
        '-ico',
        '--input_cl_obo',
        type=str,
        required=True,
        help='Input EBI Cell Ontology OBO file'
    )
    parser.add_argument(
        '-oc',
        '--output_cells',
        required=True,
        help='Output AnnData-Zarr store'
    )
    parser.add_argument(
        '-ocs',
        '--output_cell_sets',
        required=True,
        help='Output obsSets.json file'
    )
    args = parser.parse_args()
    convert_to_zarr_and_json(
        args.input_cells_h5ad,
        args.input_annotations_csv,
        args.input_cl_obo,
        args.output_cells,
        args.output_cell_sets
    )
