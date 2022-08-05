from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import pandas as pd
import numpy as np


def create_test_anndata_file(h5ad_path):
    X = np.array([
        [0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 1.6466103],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 1.8969784, 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 1.8969784, 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 2.0868707, 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 1.9842411],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 1.6075206, 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 1.6075206, 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0.],
        [0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.,
            0., 0., 2.0462136, 0., 0., 0.,
            0., 0.]
    ])

    var_arr = [
        'RP11-34P13.7', 'AP006222.2', 'RP4-669L17.10', 'RP5-857K21.4',
        'RP11-206L10.9', 'LINC00115', 'RP11-54O7.1', 'LINC02593', 'SAMD11',
        'NOC2L', 'KLHL17', 'PLEKHN1', 'PERM1', 'RP11-54O7.17', 'HES4',
        'ISG15', 'RP11-54O7.11', 'AGRN', 'RNF223', 'C1orf159'
    ]
    var_df = pd.DataFrame(
        data=[
            {
                'index': i,
                'highly_variable': True
            }
            for i in var_arr
        ]
    )
    var_df = var_df.set_index('index')

    obs_index_arr = [
        'hHP1_ACTCAATAGCAA-habib17', 'hHP1_TTCCCGTTAAAG-habib17',
        'hHP1_GTCATTGAATCA-habib17', 'hHP1_CACCTTCAATAC-habib17',
        'hHP1_ATACATGTTGTC-habib17', 'hHP1_GTTTATTAATGG-habib17',
        'hHP1_GAGTGTGAACAN-habib17', 'hHP1_TATTCCTGAGAA-habib17',
        'hHP1_ACCCAACTTCAG-habib17', 'hHP1_TATGCTCGAATT-habib17',
        'hHP1_TTCTAAACCGAC-habib17', 'hHP1_TCGGTTGCCCAT-habib17',
        'hHP1_GTCTAGGGTCTC-habib17', 'hHP1_CCTAGCATTAGN-habib17',
        'hHP1_TTAAAGAAGGAG-habib17', 'hHP1_CGCTCTAGTATC-habib17',
        'hHP1_ATCCGCAGCCGN-habib17', 'hHP1_TGGCGCAGAAGG-habib17',
        'hHP1_TATTCACCGCTT-habib17', 'hHP1_GACCATCTTATT-habib17'
    ]
    obs_celltype_arr = [
        'exCA1',
        'exCA3',
        'ASC1',
        'exCA1',
        'exCA3',
        'GABA1',
        'exCA1',
        'exCA3',
        'GABA1',
        'ODC1',
        'exCA1',
        'exCA1',
        'exDG',
        'ODC1',
        'Unclassified',
        'exDG',
        'GABA1',
        'exCA1',
        'exPFC2',
        'GABA2'
    ]
    obs_df = pd.DataFrame(
        data=[
            {'index': i, 'CellType': ct}
            for i, ct in zip(obs_index_arr, obs_celltype_arr)
        ]
    )
    obsm = {"X_umap": np.array([[0, 1] for c in obs_index_arr])}
    try:
        from anndata import AnnData
    except ModuleNotFoundError:  # pragma: no cover
        # TODO: When we don't need backward compatibility, move import back to top.
        return
    adata = AnnData(X, var=var_df, obs=obs_df, obsm=obsm)
    adata.write_h5ad(h5ad_path)


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
        '2:55001-60000',
        '3:15001-20000',
        '18:10001-15000',
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
    barcodes_df.to_csv(barcodes_path, sep='\t', index=False, header=False)
    clusters_df.to_csv(clusters_path, index=True)

    mtx = np.array([
        [0, 2, 1, 3, 0, 4, 9, 0, 0, 1, 0, 0, 0, 1, 0, 3],
        [1, 1, 3, 1, 0, 0, 0, 0, 0, 2, 2, 3, 4, 2, 4, 3],
        [0, 1, 1, 1, 1, 1, 3, 2, 1, 0, 0, 0, 1, 3, 5, 0],
        [0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 2, 1, 3, 1, 0, 0],
        [2, 3, 2, 4, 1, 2, 3, 1, 0, 0, 0, 0, 1, 0, 0, 0],
        [4, 0, 4, 0, 1, 0, 0, 0, 8, 1, 4, 3, 2, 1, 0, 2],
    ])
    coo_mtx = coo_matrix(mtx)

    mmwrite(mtx_path, coo_mtx)
