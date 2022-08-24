import argparse
import json
from anndata import read_h5ad
import pandas as pd
import numpy as np
import scipy.cluster

from constants import COLUMNS
from generate_cell_sets import (
    generate_leiden_cluster_cell_sets,
    generate_cell_type_cell_sets
)
from utils import (
    merge_cell_sets_trees
)


def clean_gexp(adata):
    gexp_arr = adata.X
    gexp_df = adata.to_df()

    # Re-scale the gene expression values between 0 and 255 (one byte ints).
    gexp_arr_min = gexp_arr.min()
    gexp_arr_max = gexp_arr.max()
    gexp_arr_range = gexp_arr_max - gexp_arr_min
    gexp_arr_ratio = 255 / gexp_arr_range
    gexp_norm_arr = (gexp_arr - gexp_arr_min) * gexp_arr_ratio

    # Perform hierarchical clustering along the genes axis.
    Z = scipy.cluster.hierarchy.linkage(gexp_norm_arr.T, method="ward")
    labels = adata.var.index.values

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = labels[leaf_index_list].tolist()

    # Create a new *ordered* gene expression dataframe.
    gexp_norm_df = pd.DataFrame(
        index=gexp_df.index.values.tolist(),
        columns=gexp_df.columns.values.tolist(),
        data=gexp_norm_arr
    )
    adata.X = gexp_norm_df.values
    
    # Sort by selecting along the gene axis of the AnnData object
    adata_sorted = adata[:, leaf_list].copy()
    adata.X = adata.X.astype(np.dtype('uint8'))
    return adata_sorted


def generate_json_files(
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
    cell_type_cell_sets = generate_cell_type_cell_sets(
        df,
        input_cl_obo_file
    )

    # Merge the Leiden Cluster and Cell Type Annotation cell sets.
    cell_sets = merge_cell_sets_trees(
        leiden_cell_sets,
        cell_type_cell_sets
    )

    with open(output_cell_sets, 'w') as f:
        json.dump(cell_sets, f, indent=1)

    # Round prediction score for labels.
    adata.obs[COLUMNS.PREDICTION_SCORE.value] = adata.obs[COLUMNS.PREDICTION_SCORE.value].apply(lambda val: '{:.2g}'.format(val))
    
    adata = clean_gexp(adata)
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
    generate_json_files(
        args.input_cells_h5ad,
        args.input_annotations_csv,
        args.input_cl_obo,
        args.output_cells,
        args.output_cell_sets
    )
