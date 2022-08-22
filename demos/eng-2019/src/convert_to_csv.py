import argparse
import json
import pandas as pd


def convert_to_csv(args):
    with open(args.input_cells) as f:
        cells_json = json.load(f)

    cell_ids = cells_json.keys()
    cells_df = pd.DataFrame(index=cell_ids, columns=["TSNE_1", "TSNE_2", "UMAP_1", "UMAP_2", "Leiden", "Kmeans", "X", "Y"])
    cells_df['TSNE_1'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['t-SNE'][0], axis='columns')
    cells_df['TSNE_2'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['t-SNE'][1], axis='columns')
    cells_df['UMAP_1'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['UMAP'][0], axis='columns')
    cells_df['UMAP_2'] = cells_df.apply(lambda row: cells_json[str(row.name)]['mappings']['UMAP'][1], axis='columns')
    cells_df['X'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][0], axis='columns')
    cells_df['Y'] = cells_df.apply(lambda row: cells_json[str(row.name)]['xy'][1], axis='columns')
    cells_df['Leiden'] = cells_df.apply(lambda row: cells_json[str(row.name)]['factors']['pleiden_clus'], axis='columns')
    cells_df['Kmeans'] = cells_df.apply(lambda row: cells_json[str(row.name)]['factors']['kmeans'], axis='columns')
    cells_df.index = cells_df.index.rename('cell_id')

    def to_diamond(x, y, r):
        return [[x, y + r], [x + r, y], [x, y - r], [x - r, y]]

    segmentations = dict()
    for cell_id, cell_obj in cells_json.items():
        segmentations[cell_id] = to_diamond(cell_obj['xy'][0], cell_obj['xy'][1], 50)

    cells_df.to_csv(args.output_cells, index=True)

    with open(args.output_segmentations, 'w') as f:
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
        '-oc',
        '--output_cells',
        type=str,
        required=True,
        help='Output cells.csv file'
    )
    parser.add_argument(
        '-os',
        '--output_segmentations',
        type=str,
        required=True,
        help='Output obsSegmentations.json file'
    )
    args = parser.parse_args()
    convert_to_csv(args)
