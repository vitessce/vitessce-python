import argparse

import scanpy as sc

def create_zarr(output_path):
    adata = sc.datasets.visium_sge(sample_id="V1_Adult_Mouse_Brain", include_hires_tiff=True)
    adata.X = adata.X.toarray()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype('uint8')
    adata.write_zarr(output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='Output Zarr store'
    )
    args = parser.parse_args()
    create_zarr(
        args.output
    )
