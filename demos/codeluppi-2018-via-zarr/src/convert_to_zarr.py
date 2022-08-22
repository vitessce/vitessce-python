import argparse
import json
import pandas as pd
import numpy as np
from anndata import AnnData


def to_uint8(arr):
    arr *= 255.0 / arr.max()
    arr = arr.astype(np.dtype('uint8'))
    return arr


def convert_to_csv(args):
    with open(args.input_cells) as f:
        cells_json = json.load(f)
    with open(args.input_molecules) as f:
        molecules_json = json.load(f)
    with open(args.input_expression_matrix) as f:
        expression_matrix = json.load(f)

    molecules_df = pd.DataFrame(index=[], columns=["X", "Y", "Gene"])
    start_at = 0
    for molecule_type, molecule_arr in molecules_json.items():
        molecule_type_ids = [str(i) for i in range(start_at, start_at + len(molecule_arr))]
        start_at += len(molecule_arr)
        molecule_type_df = pd.DataFrame(index=molecule_type_ids, columns=molecules_df.columns.values.tolist())
        molecule_type_df['Gene'] = molecule_type
        molecule_type_df['X'] = [m[0] for m in molecule_arr]
        molecule_type_df['Y'] = [m[1] for m in molecule_arr]

        molecules_df = pd.concat([molecules_df, molecule_type_df])
    molecules_df.index = molecules_df.index.rename('molecule_id')

    cell_ids = cells_json.keys()
    cells_df = pd.DataFrame(index=cell_ids, columns=["TSNE_1", "TSNE_2", "PCA_1", "PCA_2", "Cluster", "Subcluster", "X", "Y"])
    cells_df['TSNE_1'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['t-SNE'][0], axis='columns')
    cells_df['TSNE_2'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['t-SNE'][1], axis='columns')
    cells_df['PCA_1'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['PCA'][0], axis='columns')
    cells_df['PCA_2'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['PCA'][1], axis='columns')
    cells_df['X'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][0], axis='columns')
    cells_df['Y'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][1], axis='columns')
    cells_df['Cluster'] = cells_df.apply(lambda row: cells_json[str(row.name)]['factors']['cluster'], axis='columns')
    cells_df['Subcluster'] = cells_df.apply(lambda row: cells_json[str(row.name)]['factors']['subcluster'], axis='columns')
    cells_df.index = cells_df.index.rename('cell_id')

    matrix_df = pd.DataFrame(index=expression_matrix['cols'], columns=expression_matrix['rows'])
    for i, gene in enumerate(matrix_df.columns.values.tolist()):
        matrix_df[gene] = expression_matrix['matrix'][i]
    matrix_df.index = matrix_df.index.rename('cell_id')

    # Cells zarr
    c_obs_df = cells_df[['Cluster', 'Subcluster']]
    c_var_df = matrix_df.T[[]]

    c_adata = AnnData(X=matrix_df.values, obs=c_obs_df, var=c_var_df)
    c_adata.obsm['X_tsne'] = cells_df[['TSNE_1', 'TSNE_2']].values.astype(np.dtype('<f4'))
    c_adata.obsm['X_pca'] = cells_df[['PCA_1', 'PCA_2']].values.astype(np.dtype('<f4'))
    c_adata.obsm['X_spatial'] = cells_df[['X', 'Y']].values.astype(np.dtype('<u4'))
    c_adata.X = to_uint8(c_adata.X)

    num_cells = len(cells_json.keys())
    c_adata.obsm['X_segmentations'] = np.zeros((num_cells, len(list(cells_json.values())[0]['poly']), 2))
    for i, cell_id in enumerate(c_adata.obs.index.values.tolist()):
        c_adata.obsm['X_segmentations'][i, :, :] = cells_json[cell_id]['poly']
    c_adata.obsm['X_segmentations'] = c_adata.obsm['X_segmentations'].astype(np.dtype('<u4'))

    # Molecules zarr
    m_obs_df = molecules_df[['Gene']]
    m_var_df = pd.DataFrame(index=[], columns=[], data=[])

    m_adata = AnnData(X=None, obs=m_obs_df, var=m_var_df)
    m_adata.obsm['X_spatial'] = molecules_df[['X', 'Y']].values.astype(np.dtype('<u4'))

    # Write outputs
    c_adata.write_zarr(args.output_cells)
    m_adata.write_zarr(args.output_molecules)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-ic',
        '--input_cells',
        type=str,
        required=True,
        help='Input cells.json file'
    )
    parser.add_argument(
        '-im',
        '--input_molecules',
        type=str,
        required=True,
        help='Input molecules.json file'
    )
    parser.add_argument(
        '-iem',
        '--input_expression_matrix',
        type=str,
        required=True,
        help='Input clusters.json file'
    )
    parser.add_argument(
        '-oc',
        '--output_cells',
        type=str,
        required=True,
        help='Output zarr store for cells'
    )
    parser.add_argument(
        '-om',
        '--output_molecules',
        type=str,
        required=True,
        help='Output zarr store for molecules'
    )
    args = parser.parse_args()
    convert_to_csv(args)
