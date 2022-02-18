import argparse
import json
from anndata import AnnData
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
import pickle

def convert_to_zarr(input_molecules, input_polygon, input_gene_colors, input_gene_array, output_path):
    molecules_table = pq.read_table(input_molecules)
    df = molecules_table.to_pandas()
    df.index = df.index.rename("molecule_id")

    # TODO: remove
    df = df.head(500000)
    
    with open(input_gene_colors, "rb") as f:
        gene_colors = pickle.load(f)

    df["r"] = df.apply(lambda row: gene_colors[row["decoded_genes"]][0], axis='columns')
    df["g"] = df.apply(lambda row: gene_colors[row["decoded_genes"]][1], axis='columns')
    df["b"] = df.apply(lambda row: gene_colors[row["decoded_genes"]][2], axis='columns')

    rgb_df = df[["r", "g", "b"]]
    rgb_df *= 255

    rgb_df = rgb_df.apply(lambda col: np.floor(col), axis='rows')

    obsm = {
        "spatial": df[["r_px_global_stitched", "c_px_global_stitched"]].values.astype('int16'),
        "rgb": rgb_df.values.astype('uint8'),
    }

    X = np.zeros((df.shape[0], 0))
    var_df = pd.DataFrame(data=[])
    obs_df = pd.DataFrame(index=df.index.values.tolist(), data=df[["decoded_genes"]].values, columns=["gene_id"])
    #obs_df["cell_id"] = obs_df["cell_id"].astype(np.uint16)

    gene_names = list(df["decoded_genes"].unique())
    obs_df["gene_index"] = obs_df["gene_id"].apply(lambda gene_name: gene_names.index(gene_name))
    obs_df["gene_index"] = obs_df["gene_index"].astype(np.uint16)

    mol_adata = AnnData(X=None, obs=obs_df, obsm=obsm, var=None)
    mol_adata.write_zarr(output_path)

    print(mol_adata)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-im',
        '--input_molecules',
        type=str,
        required=True,
        help='Input parquet file'
    )
    parser.add_argument(
        '-ip',
        '--input_polygon',
        type=str,
        required=True,
        help='Input csv file'
    )
    parser.add_argument(
        '-igc',
        '--input_gene_colors',
        type=str,
        required=True,
        help='Input pkl file'
    )
    parser.add_argument(
        '-iga',
        '--input_gene_array',
        type=str,
        required=True,
        help='Input pkl file'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='Output Zarr store'
    )
    args = parser.parse_args()
    convert_to_zarr(
        args.input_molecules,
        args.input_polygon,
        args.input_gene_colors,
        args.input_gene_array,
        args.output
    )
