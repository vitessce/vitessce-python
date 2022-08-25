import unittest
from anndata import AnnData
import pandas as pd
import numpy as np

from vitessce.data_utils import (
    optimize_arr,
    optimize_adata,
    sort_var_axis,
    to_uint8,
)


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
        assert norm_X.tolist() == [[255, 200, 250], [250, 50, 255], [104, 50, 100], [200, 70, 205]]

    def test_to_uint8_gene_norm(self):
        adata = self.adata
        norm_X = to_uint8(adata.X, norm_along="var")
        assert norm_X.tolist() == [[255, 255, 246], [246, 0, 254], [0, 0, 0], [161, 33, 172]]

    def test_to_uint8_cell_norm(self):
        adata = self.adata
        norm_X = to_uint8(adata.X, norm_along="obs")
        assert norm_X.tolist() == [[255, 0, 231], [248, 0, 255], [255, 0, 231], [245, 0, 255]]
