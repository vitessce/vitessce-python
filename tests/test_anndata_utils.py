import unittest

from anndata import AnnData
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.io import mmread
import zarr

from vitessce.data_utils import (
    optimize_arr,
    optimize_adata,
    sort_var_axis,
    to_uint8,
    adata_to_multivec_zarr,
)


data_path = Path('tests/data')


class TestAnnDataUtils(unittest.TestCase):

    def setUp(self):
        X = np.array([
            [5.1, 4.0, 5.0],
            [5.0, 1.0, 5.1],
            [2.1, 1.0, 2.0],
            [4.0, 1.4, 4.1],
        ])
        obs = pd.DataFrame(data=[
            {"cell_id": "cell_1", "leiden": "Cluster 2", "bad": "_"},
            {"cell_id": "cell_2", "leiden": "Cluster 2", "bad": "+"},
            {"cell_id": "cell_3", "leiden": "Cluster 1", "bad": "%"},
            {"cell_id": "cell_4", "leiden": "Cluster 5", "bad": "$"}
        ]).set_index("cell_id")
        var = pd.DataFrame(data=[
            {"gene_id": "SOX10"},
            {"gene_id": "LAMP5"},
            {"gene_id": "RORB"}
        ]).set_index("gene_id")
        obsm = {
            "X_umap": np.array([
                [0.5, 1.0],
                [1.0, 2.0],
                [3.1, 2.2],
                [4.0, 1.0]
            ]),
            "X_ints": np.array([
                [1.0, 2.0],
                [1.0, 2.0],
                [3.0, 2.0],
                [4.0, 1.0]
            ])
        }
        self.adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)

    def test_optimize_arr(self):
        int_arr = self.adata.obsm["X_ints"]
        assert optimize_arr(int_arr).dtype.name == 'uint8'

        float_arr = self.adata.obsm["X_umap"]
        assert optimize_arr(float_arr).dtype.name == 'float32'

    def test_optimize_adata(self):
        adata = self.adata
        adata_optim = optimize_adata(adata, obs_cols=['leiden'], obsm_keys=['X_umap'])
        assert list(adata_optim.obsm.keys()) == ['X_umap']
        assert adata_optim.obsm['X_umap'].dtype.name == 'float32'
        assert adata_optim.obs.columns.values.tolist() == ['leiden']

    def test_sort_var_axis(self):
        adata = self.adata
        leaf_list = sort_var_axis(adata.X, adata.var.index.values)
        adata_sorted = adata[:, leaf_list].copy()
        assert adata.var.index.values.tolist() == ['SOX10', 'LAMP5', 'RORB']
        assert adata_sorted.var.index.values.tolist() == ['LAMP5', 'SOX10', 'RORB']

    def test_to_uint8(self):
        adata = self.adata
        norm_X = to_uint8(adata.X)
        assert norm_X.tolist() == [[5, 4, 5], [5, 1, 5], [2, 1, 2], [4, 1, 4]]

    def test_to_uint8_global_norm(self):
        adata = self.adata
        norm_X = to_uint8(adata.X, norm_along="global")
        np.testing.assert_almost_equal(norm_X.tolist(), [[255, 200, 250], [250, 50, 255], [104, 50, 100], [200, 70, 205]], decimal=0)

    def test_to_uint8_gene_norm(self):
        adata = self.adata
        norm_X = to_uint8(adata.X, norm_along="var")
        np.testing.assert_almost_equal(norm_X.tolist(), [[255, 255, 246], [246, 0, 254], [0, 0, 0], [161, 33, 172]], decimal=0)

    def test_to_uint8_cell_norm(self):
        adata = self.adata
        norm_X = to_uint8(adata.X, norm_along="obs")
        assert norm_X.tolist() == [[255, 0, 231], [248, 0, 255], [255, 0, 231], [245, 0, 255]]

    def test_multivec_zarr(self):
        mtx = mmread(data_path / 'test.snap.mtx').toarray()
        bins_df = pd.read_csv(
            data_path / 'test.snap.bins.txt', header=None, names=["interval"])
        clusters_df = pd.read_csv(
            data_path / 'test.snap.clusters.csv', index_col=0)

        zarr_filepath = data_path / 'test_out.snap.multivec.zarr'

        # The genome assembly is GRCh38 but the chromosome names in the bin names do not start with the "chr" prefix.
        # This is incompatible with the chromosome names from `negspy`, so we need to append the prefix.
        bins_df["interval"] = bins_df["interval"].apply(lambda x: "chr" + x)

        obs = clusters_df[["cluster"]]
        obs["cluster"] = obs["cluster"].astype(str)
        obsm = {"X_umap": clusters_df[["umap.1", "umap.2"]].values}
        adata = AnnData(X=mtx, obs=obs, var=bins_df, obsm=obsm)

        # Sort cluster IDs
        cluster_ids = obs["cluster"].unique().tolist()
        cluster_ids.sort(key=int)
        # Save genomic profiles to multivec-zarr format.
        adata_to_multivec_zarr(adata, zarr_filepath, obs_set_col="cluster", obs_set_name="Cluster", obs_set_vals=cluster_ids)

        z = zarr.open(zarr_filepath, mode='r')

        self.assertEqual(z['chromosomes/chr1/5000'].shape, (4, 49792))
        self.assertEqual(z['chromosomes/chr1/5000'][:, 0].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][:, 1].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][:, 2].sum(), 7)
        self.assertEqual(z['chromosomes/chr1/5000'][:, 3].sum(), 7)
        self.assertEqual(z['chromosomes/chr1/5000'][:, 4].sum(), 0)
        self.assertEqual(z['chromosomes/chr1/5000'][0, :].sum(), 17)
        self.assertEqual(z['chromosomes/chr1/10000'][0, :].sum(), 17)
        self.assertEqual(z['chromosomes/chr1/5000'][:, 2].sum(), 7)
        self.assertEqual(z['chromosomes/chr2/5000'][:, 0].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/5000'][:, 1].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/10000'][:, 0].sum(), 0)
        self.assertEqual(z['chromosomes/chr2/5000'][:, 2].sum(), 4)
        self.assertEqual(z['chromosomes/chr2/5000'][:, 3].sum(), 9)
        self.assertEqual(z['chromosomes/chr2/10000'][:, 1].sum(), 13)
        self.assertEqual(z['chromosomes/chr3/5000'][:, 3].sum(), 9)
        self.assertEqual(z['chromosomes/chr3/5000'][:].sum(), 9)
        self.assertEqual(z['chromosomes/chr18/5000'][:].sum(), 8)
