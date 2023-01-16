import argparse
import scanpy as sc
import numpy as np
import scipy.cluster
from vitessce.data_utils import (
    to_diamond,
    rgb_img_to_ome_zarr,
    optimize_adata,
)


def create_zarr(output_adata, output_img):
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
    # Determine the top 300 highly variable genes.
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=300)

    # Dimensionality reduction and clustering
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters")

    # Hierarchical clustering of genes for optimal gene ordering
    X_hvg_arr = adata[:, adata.var['highly_variable']].X.toarray()
    X_hvg_index = adata[:, adata.var['highly_variable']].var.copy().index

    Z = scipy.cluster.hierarchy.linkage(X_hvg_arr.T, method="average", optimal_ordering=True)

    # Get the hierarchy-based ordering of genes.
    num_cells = adata.obs.shape[0]
    highly_var_index_ordering = scipy.cluster.hierarchy.leaves_list(Z)
    highly_var_genes = X_hvg_index.values[highly_var_index_ordering].tolist()

    all_genes = adata.var.index.values.tolist()
    not_var_genes = adata.var.loc[~adata.var['highly_variable']].index.values.tolist()

    def get_orig_index(gene_id):
        return all_genes.index(gene_id)

    var_index_ordering = list(map(get_orig_index, highly_var_genes)) + list(map(get_orig_index, not_var_genes))

    # Create a new *ordered* gene expression dataframe.
    adata = adata[:, var_index_ordering].copy()
    adata.obsm["X_hvg"] = adata[:, adata.var['highly_variable']].X.copy()

    # Unclear what the exact scale factor is required to align
    # the spots to the image. Through trial and error / manual binary search
    # of values I arrived at:
    scale_factor = 1 / 5.87
    adata.obsm['spatial'] = (adata.obsm['spatial'] * scale_factor)

    adata.obsm['segmentations'] = np.zeros((num_cells, 4, 2))
    radius = 10
    for i in range(num_cells):
        adata.obsm['segmentations'][i, :, :] = to_diamond(adata.obsm['spatial'][i, 0], adata.obsm['spatial'][i, 1], radius)

    # Write img_arr to OME-Zarr.
    # Need to convert images from interleaved to non-interleaved (color axis should be first).
    img_hires = adata.uns['spatial']['V1_Human_Lymph_Node']['images']['hires']
    img_arr = np.transpose(img_hires, (2, 0, 1))
    # Convert values from [0, 1] to [0, 255].
    img_arr *= 255.0

    rgb_img_to_ome_zarr(img_arr, output_img, axes="cyx", chunks=(1, 256, 256), img_name="H & E Image")
    adata = optimize_adata(
        adata,
        obs_cols=["clusters"],
        var_cols=["highly_variable"],
        obsm_keys=["X_hvg", "spatial", "segmentations", "X_umap", "X_pca"],
        optimize_X=True,
        # Vitessce plays nicely with dense matrices saved with chunking
        # and this one is small enough that dense is not a huge overhead.
        to_dense_X=True,
    )
    adata.write_zarr(output_adata, chunks=[adata.shape[0], 10])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-oa',
        '--output_adata',
        type=str,
        required=True,
        help='Output AnnData-Zarr store'
    )
    parser.add_argument(
        '-oi',
        '--output_img',
        type=str,
        required=True,
        help='Output OME-Zarr store'
    )
    args = parser.parse_args()
    create_zarr(
        args.output_adata,
        args.output_img
    )
