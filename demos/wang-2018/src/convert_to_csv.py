import argparse
import json
import pandas as pd
from pathlib import Path

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
    cells_df = pd.DataFrame(index=cell_ids, columns=["X", "Y"])
    cells_df['X'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][0], axis='columns')
    cells_df['Y'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][1], axis='columns')
    cells_df.index = cells_df.index.rename('cell_id')

    segmentations = dict()
    for cell_id, cell_obj in cells_json.items():
        segmentations[cell_id] = cell_obj['poly']

    matrix_df = pd.DataFrame(index=cell_ids, columns=list(expression_matrix.keys()))
    matrix_df.index = matrix_df.index.astype(str)
    for gene_id in matrix_df.columns.values.tolist():
        for cell_id in matrix_df.index.values.tolist():
            try:
                matrix_df.at[cell_id, gene_id] = expression_matrix[gene_id]['cells'][cell_id] / expression_matrix[gene_id]['max']
            except KeyError:
                matrix_df.at[cell_id, gene_id] = 0
    matrix_df.index = matrix_df.index.rename('cell_id')

    # Pandas does not make intermediate paths
    output_dir = Path(args.output_cells).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    cells_df.to_csv(args.output_cells, index=True)
    matrix_df.to_csv(args.output_cells_matrix, index=True)
    molecules_df.to_csv(args.output_molecules, index=True)

    with open(args.output_cells_segmentations, 'w') as f:
        json.dump(segmentations, f)


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
        help='Output cells.csv file'
    )
    parser.add_argument(
        '-ocs',
        '--output_cells_segmentations',
        type=str,
        required=True,
        help='Output obsSegmentations.json file'
    )
    parser.add_argument(
        '-ocm',
        '--output_cells_matrix',
        type=str,
        required=True,
        help='Output obsFeatureMatrix.csv file'
    )
    parser.add_argument(
        '-om',
        '--output_molecules',
        type=str,
        required=True,
        help='Output molecules.csv file'
    )
    args = parser.parse_args()
    convert_to_csv(args)
