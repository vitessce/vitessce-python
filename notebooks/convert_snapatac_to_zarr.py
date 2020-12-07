from vitessce import SnapWrapper
from os.path import join
from scipy.io import mmread
import pandas as pd
import numpy as np

if __name__ == "__main__":
    mtx = mmread(join('data', 'snapatac', 'filtered_cell_by_bin.mtx'))
    barcodes_df = pd.read_csv(join('data', 'snapatac', 'barcodes.txt'), header=None)
    bins_df = pd.read_csv(join('data', 'snapatac', 'bins.txt'), header=None)
    clusters_df = pd.read_csv(join('data', 'snapatac', 'umap_coords_clusters.csv'), index_col=0)

    zarr_filepath = join('data', 'snapatac', 'out.snap.multivec.zarr')

    w = SnapWrapper(mtx, barcodes_df, bins_df, clusters_df)
    w.create_genomic_multivec_zarr(zarr_filepath)

    print("Done")
