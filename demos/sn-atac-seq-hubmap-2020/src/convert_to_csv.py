import argparse
import json
import pandas as pd


def convert_to_csv(args):
    with open(args.input_cells) as f:
        cells_json = json.load(f)
    with open(args.input_cell_sets) as f:
        cell_sets_json = json.load(f)

    cell_ids = cells_json.keys()
    cells_df = pd.DataFrame(index=cell_ids, columns=["UMAP_1", "UMAP_2", "cluster"])
    cells_df['UMAP_1'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['UMAP'][0], axis='columns')
    cells_df['UMAP_2'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['UMAP'][1], axis='columns')
    cells_df.index = cells_df.index.rename('cell_id')

    for cluster in cell_sets_json['tree'][0]['children']:
        cluster_name = cluster['name']
        cluster_cell_ids = cluster['set']

        for cell_id in cluster_cell_ids:
            cells_df.at[cell_id, 'cluster'] = cluster_name

    cells_df['cluster'] = cells_df['cluster'].astype(str)

    cells_df.to_csv(args.output_cells, index=True)


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
        '-ics',
        '--input_cell_sets',
        type=str,
        required=True,
        help='Input molecules.json file'
    )
    parser.add_argument(
        '-oc',
        '--output_cells',
        type=str,
        required=True,
        help='Output cells.csv file'
    )
    args = parser.parse_args()
    convert_to_csv(args)
