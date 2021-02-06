import argparse

from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.cluster

def create_zarr(output_path):
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node", include_hires_tiff=True)

    # Reference: https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html
    # Calculate QC metrics
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Perform basic filtering
    sc.pp.filter_cells(adata, min_counts=5000)
    sc.pp.filter_cells(adata, max_counts=35000)
    adata = adata[adata.obs["pct_counts_mt"] < 20]
    sc.pp.filter_genes(adata, min_cells=10)



    # Perform normalization
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    # Determine the top 2000 highly variable genes.
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=300)


    # Hierarchical clustering of genes for optimal gene ordering
    X_hvg_arr = adata[:, adata.var['highly_variable']].X.copy().toarray()
    X_hvg_index = adata[:, adata.var['highly_variable']].var.copy().index

    Z = scipy.cluster.hierarchy.linkage(X_hvg_arr.T, method="average", optimal_ordering=True)
    labels = X_hvg_index.values

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = labels[leaf_index_list].tolist()

    # Create a new *ordered* gene expression dataframe.
    X_hvg_df = pd.DataFrame(
        index=adata.obs.index.values.tolist(),
        columns=X_hvg_index.values.tolist(),
        data=X_hvg_arr
    )
    X_hvg_ordered_df = X_hvg_df[leaf_list]

    highly_variable_ordered = leaf_list
    not_variable = set(adata.var.index.values.tolist()) - set(highly_variable_ordered)
    highly_variable_ordered_plus_not_variable = highly_variable_ordered + list(not_variable)

    X_df = pd.DataFrame(
        index=adata.obs.index.values.tolist(),
        columns=adata.var.index.values.tolist(),
        data=adata.X.toarray()
    )
    X_ordered_df = X_df[highly_variable_ordered_plus_not_variable]

    var_df = adata.var
    var_ordered_df = var_df.loc[highly_variable_ordered_plus_not_variable]

    adata_uns_spatial = adata.uns['spatial']

    # Create a new AnnData object with *ordered* genes.
    adata = AnnData(X=X_ordered_df.values, var=var_ordered_df.copy(), obs=adata.obs.copy())
    adata.obsm['X_hvg'] = X_hvg_ordered_df.values
    adata.uns['spatial'] = adata_uns_spatial

    # Dimensionality reduction and clustering
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters")

    adata.obsm['xy'] = np.stack((adata.obs['array_col'].values, adata.obs['array_row'].values), axis=-1).astype('uint8')

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
