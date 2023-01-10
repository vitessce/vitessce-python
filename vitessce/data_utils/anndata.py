import numpy as np
from anndata import AnnData
import scipy.cluster
from scipy.sparse import issparse


# Try to cast an array to a dtype that takes up less space.
def cast_arr(arr):
    orig_sum = np.sum(arr)
    orig_max = np.max(arr)
    orig_min = np.min(arr)

    # Try casting float to int for better downstream compression.
    if arr.dtype.kind == 'f':
        cast_arr = arr.astype(f'<i{arr.dtype.itemsize}')
        cast_sum = np.sum(cast_arr)
        if np.abs(orig_sum - cast_sum) < 1e-2:
            arr = cast_arr

    # Try casting signed int to unsigned int.
    if arr.dtype.kind == 'i':
        # If signed int dtype but no values are less than zero,
        # convert to unsigned int dtype.
        if arr.min() >= 0:
            arr = arr.astype(f'<u{arr.dtype.itemsize}')

    # Try casting to a smaller itemsize.
    if arr.dtype.kind == 'u' or arr.dtype.kind == 'i' or arr.dtype.kind == 'f':
        next_itemsizes = [4] if arr.dtype.kind == 'f' else [4, 2, 1]
        for next_itemsize in next_itemsizes:
            if arr.dtype.itemsize > next_itemsize:
                next_dtype = np.dtype(f'<{arr.dtype.kind}{next_itemsize}')
                next_dtype_info = np.iinfo(next_dtype) if arr.dtype.kind == 'u' or arr.dtype.kind == 'i' else np.finfo(next_dtype)
                if next_dtype_info.min <= orig_min and next_dtype_info.max >= orig_max:
                    arr = arr.astype(next_dtype)
                elif arr.dtype.itemsize == 8 and (arr.dtype.kind == 'u' or arr.dtype.kind == 'i'):
                    print(f"WARNING: Not casting array with dtype {arr.dtype.name}, but Zarr.js suggests avoiding int64 and uint64")

    # Check for float16 usage.
    if arr.dtype.kind == 'f' and arr.dtype.itemsize == 2:
        # Zarr.js does not have a Float16Array type
        arr = arr.astype('<f4')

    return arr


def optimize_arr(arr):
    arr = to_dense(cast_arr(to_memory(arr)))
    if isinstance(arr, np.matrix):
        arr = np.array(arr)
    return arr


# Given an anndata object, optimize for usage with Vitessce
def optimize_adata(adata, obs_cols=None, obsm_keys=None, var_cols=None, varm_keys=None, layer_keys=None, preserve_X=False, remove_X=False):
    if not remove_X:
        if not preserve_X and adata.X is not None:
            new_X = optimize_arr(adata.X)
        else:
            new_X = adata.X
    else:
        new_X = None
    # Only keep the subset of columns required.
    if obs_cols is not None:
        new_obs = adata.obs[obs_cols]
    else:
        new_obs = adata.obs
    if var_cols is not None:
        new_var = adata.var[var_cols]
    else:
        new_var = adata.var
    # Only keep the subset of obsm and varm items required.
    if obsm_keys is None:
        obsm_keys = adata.obsm.keys()
    new_obsm = {
        obsm_key: optimize_arr(adata.obsm[obsm_key])
        for obsm_key
        in obsm_keys
    }
    if varm_keys is None:
        varm_keys = adata.varm.keys()
    new_varm = {
        varm_key: optimize_arr(adata.varm[varm_key])
        for varm_key
        in varm_keys
    }
    if layer_keys is None:
        layer_keys = adata.layers.keys()
    new_layers = {
        layer_key: optimize_arr(adata.layers[layer_key])
        for layer_key
        in layer_keys
    }
    adata = AnnData(X=new_X, obs=new_obs, var=new_var, obsm=new_obsm, varm=new_varm, layers=new_layers)
    return adata


def to_memory(arr):
    if type(arr).__name__ == "SparseDataset":
        # This comes from using a backed AnnData object (e.g. read_h5ad(path, backed="r+")).
        # Reference: https://github.com/scverse/anndata/blob/286bc7f207863964cb861f17c96ab24fe0cf72ac/anndata/_core/sparse_dataset.py#L230
        return arr.to_memory()
    return arr


# Convert a sparse array to dense.
def to_dense(arr):
    if issparse(arr):
        return arr.todense()
    return arr


# Convert an array to uint8 dtype.
def to_uint8(arr, norm_along=None):
    # Re-scale the gene expression values between 0 and 255 (one byte ints).
    if norm_along is None:
        norm_arr = arr
    elif norm_along == "global":
        arr *= 255.0 / arr.max()
        norm_arr = arr
    elif norm_along == "var":
        # Normalize along gene axis
        arr = to_dense(arr)
        num_cells = arr.shape[0]
        min_along_genes = arr.min(axis=0)
        max_along_genes = arr.max(axis=0)
        range_per_gene = max_along_genes - min_along_genes
        ratio_per_gene = 255.0 / range_per_gene

        norm_arr = np.multiply(
            (arr - np.tile(min_along_genes, (num_cells, 1))),
            np.tile(ratio_per_gene, (num_cells, 1))
        )
    elif norm_along == "obs":
        # Normalize along cell axis
        arr = to_dense(arr)
        num_genes = arr.shape[1]
        min_along_cells = arr.min(axis=1)
        max_along_cells = arr.max(axis=1)
        range_per_cell = max_along_cells - min_along_cells
        ratio_per_cell = 255.0 / range_per_cell

        norm_arr = np.multiply(
            (arr.T - np.tile(min_along_cells, (num_genes, 1))),
            np.tile(ratio_per_cell, (num_genes, 1))
        ).T
    else:
        raise ValueError("to_uint8 received unknown norm_along value")
    return norm_arr.astype('u1')


# Sort the var index after hierarchical clustering.
def sort_var_axis(adata_X, orig_var_index, full_var_index=None):
    gexp_arr = to_dense(adata_X)

    # Perform hierarchical clustering along the genes axis.
    Z = scipy.cluster.hierarchy.linkage(gexp_arr.T, method="ward")

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = orig_var_index[leaf_index_list].tolist()

    if full_var_index:
        leaf_list = leaf_list + list(set(full_var_index).difference(set(leaf_list)))

    return leaf_list


# Convert an (x, y) coordinate to a polygon (diamond) with a given radius.
def to_diamond(x, y, r):
    return np.array([[x, y + r], [x + r, y], [x, y - r], [x - r, y]])
