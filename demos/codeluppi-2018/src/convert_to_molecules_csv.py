import argparse
from anndata import read_zarr


def convert_to_csv(args):
    adata = read_zarr(args.input)

    molecules_df = adata.obs
    molecules_df['X'] = adata.obsm['X_spatial'][:, 0]
    molecules_df['Y'] = adata.obsm['X_spatial'][:, 1]
    molecules_df.index = molecules_df.index.rename("molecule_id")
    molecules_df.to_csv(args.output)


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
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output molecules.csv file"
    )
    args = parser.parse_args()
    convert_to_csv(args)
