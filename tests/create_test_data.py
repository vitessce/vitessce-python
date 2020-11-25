from os.path import join

import zarr
from anndata import read_h5ad
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import pandas as pd
import numpy as np

def create_test_anndata_file(adata_path):
    pass

def create_test_loom_file(loom_path):
    pass

def create_test_omezarr_store(zarr_path):
    pass

def create_test_ometiff_file(tiff_path):
    pass

def create_test_snaptools_files(mtx_path, bins_path, barcodes_path, clusters_path):
    bins_arr = [
        '1:10001-15000',
        '1:15001-20000',
        '1:65001-70000',
        '1:80001-85000',
        '1:105001-110000',
        '1:115001-120000',
        '1:270001-275000',
        '2:10001-15000',
        '2:15001-20000',
        '2:20001-25000',
        '2:25001-30000',
        '2:30001-35000',
        '2:35001-40000',
    ]
    bins_df = pd.DataFrame(
        data=[
            {'bin': bin_str}
            for bin_str in bins_arr
        ]
    )

    barcodes_arr = [
        'AAACATCGAACGCTTAACGTATCA',
        'AAACATCGACACGACCACACAGAA',
        'AAACATCGAGTACAAGACAGCAGA',
        'AAACATCGATCCTGTAAACCGAGA',
        'AAACATCGCACCTTACACACAGAA',
        'AAACATCGCACTTCGAAACCGAGA',
    ]
    barcodes_df = pd.DataFrame(
        data=[
            {'barcode': barcode_str}
            for barcode_str in barcodes_arr
        ]
    )

    clusters_df = pd.DataFrame(
        data=[
            {
                'barcode': 'AAACATCGAACGCTTAACGTATCA',
                'umap.1': 0.45,
                'umap.2': 1.69,
                'cluster': '4',
            },
            {
                'barcode': 'AAACATCGACACGACCACACAGAA',
                'umap.1': -1.27,
                'umap.2': -1.16,
                'cluster': '4',
            },
            {
                'barcode': 'AAACATCGAGTACAAGACAGCAGA',
                'umap.1': 4.43,
                'umap.2': 1.64,
                'cluster': '10',
            },
            {
                'barcode': 'AAACATCGATCCTGTAAACCGAGA',
                'umap.1': -0.84,
                'umap.2': 1.57,
                'cluster': '3',
            },
            {
                'barcode': 'AAACATCGCACCTTACACACAGAA',
                'umap.1': 0.54,
                'umap.2': 0.11,
                'cluster': '2',
            },
            {
                'barcode': 'AAACATCGCACTTCGAAACCGAGA',
                'umap.1': 1.24,
                'umap.2': 1.43,
                'cluster': '3',
            },
        ]
    )
    clusters_df = clusters_df.set_index('barcode')

    bins_df.to_csv(bins_path, sep='\t', index=False, header=False)
    barcodes_df.to_csv(bins_path, sep='\t', index=False, header=False)
    clusters_df.to_csv(clusters_path, index=True)

    mtx = [
        [0, 2, 1, 3, 0, 4, 9, 0, 0, 1, 0, 0, 0],
        [1, 1, 3, 1, 0, 0, 0, 0, 0, 2, 2, 3, 4],
        [0, 1, 1, 1, 1, 1, 3, 2, 1, 0, 0, 0, 1],
        [0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 2, 1, 3],
        [2, 3, 2, 4, 1, 2, 3, 1, 0, 0, 0, 0, 1],
        [4, 0, 4, 0, 1, 0, 0, 0, 8, 1, 4, 3, 2],
    ]
    coo_mtx = coo_matrix(mtx)

    mmwrite(mtx_path, coo_mtx)