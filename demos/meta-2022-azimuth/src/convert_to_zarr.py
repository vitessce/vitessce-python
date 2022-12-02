import argparse
from anndata import read_h5ad
from scipy import sparse


def convert_h5ad_to_zarr(input_path, output_path):
    adata = read_h5ad(input_path)
    # Vitessce plays nicely with csc matrices
    if isinstance(adata, sparse.spmatrix):
        adata.X = adata.X.tocsc()
    adata.write_zarr(output_path, chunks=[adata.shape[0], 10])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Input H5AD file'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='Output Zarr store'
    )
    args = parser.parse_args()
    convert_h5ad_to_zarr(
        args.input,
        args.output
    )
