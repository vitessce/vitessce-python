import numpy as np
from anndata import AnnData
import scipy.cluster
from scipy.sparse import issparse

VAR_CHUNK_SIZE = 10


def generate_h5ad_ref_spec(h5_url, omit_url=True):
    from kerchunk.hdf import SingleHdf5ToZarr
    h5chunks = SingleHdf5ToZarr(h5_url, inline_threshold=300)
    h5dict = h5chunks.translate()
    if omit_url:
        for key, val in h5dict['refs'].items():
            if isinstance(val, list):
                h5dict['refs'][key] = [None, *val[1:]]
    return h5dict


def cast_arr(arr):
    """
    Try to cast an array to a dtype that takes up less space.

    :param arr: The array to cast.
    :type arr: np.array

    :returns: The new array.
    :rtype: np.array
    """
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
    """
    Try to cast an array to a dtype that takes up less space, and convert to dense.

    :param arr: The array to cast and convert.
    :type arr: np.array

    :returns: The new array.
    :rtype: np.array
    """
    arr = to_dense(cast_arr(to_memory(arr)))
    if isinstance(arr, np.matrix):
        arr = np.array(arr)
    return arr


def optimize_adata(adata, obs_cols=None, obsm_keys=None, var_cols=None, varm_keys=None, layer_keys=None, remove_X=False, optimize_X=False, to_dense_X=False, to_sparse_X=False):
    """
    Given an AnnData object, optimize for usage with Vitessce and return a new object.

    :param adata: The AnnData object to optimize.
    :type adata: anndata.AnnData
    :param obs_cols: Columns of adata.obs to optimize. Columns not specified will not be included in the returned object.
    :type obs_cols: list[str] or None
    :param var_cols: Columns of adata.var to optimize. Columns not specified will not be included in the returned object.
    :type var_cols: list[str] or None
    :param obsm_keys: Arrays within adata.obsm to optimize. Keys not specified will not be included in the returned object.
    :type obsm_keys: list[str] or None
    :param varm_keys: Arrays within adata.varm to optimize. Keys not specified will not be included in the returned object.
    :type varm_keys: list[str] or None
    :param layer_keys: Arrays within adata.layers to optimize. Keys not specified will not be included in the returned object.
    :type layer_keys: list[str] or None
    :param bool remove_X: Should the returned object have its X matrix set to None? By default, False.
    :param bool optimize_X: Should the returned object run optimize_arr on adata.X? By default, False.
    :param bool to_dense_X: Should adata.X be cast to a dense array in the returned object? By default, False.
    :param bool to_sparse_X: Should adata.X be cast to a sparse array in the returned object? By default, False.

    :returns: The new AnnData object.
    :rtype: anndata.AnnData
    """
    if not remove_X:
        if adata.X is not None:
            if optimize_X:
                new_X = optimize_arr(adata.X)
            else:
                new_X = adata.X
            # Convert to dense or sparse if asked
            if to_dense_X:
                new_X = to_dense(new_X)
            if to_sparse_X:
                new_X = new_X.tocsc()
            if to_dense_X and to_sparse_X:
                raise ValueError(
                    "Did not expect both to_dense_X and to_sparse_X to be True")
            # If the user has not asked for their matrix to be converted to dense or sparse,
            # we still want to ensure that, if the original matrix was sparse, it
            # uses the CSC sparse format.
            if not to_dense_X and not to_sparse_X:
                # In the future, we can use sparse matrices with equal performance:
                # https://github.com/theislab/anndata/issues/524
                # Vitessce can use csc matrices somewhat efficiently.
                if issparse(new_X):
                    new_X = new_X.tocsc()
        else:
            new_X = None
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
    """
    Try to load a backed AnnData array into memory.

    :param arr: The array to load.
    :type arr: np.array

    :returns: The loaded array.
    :rtype: np.array
    """
    if type(arr).__name__ == "SparseDataset":
        # This comes from using a backed AnnData object (e.g. read_h5ad(path, backed="r+")).
        # Reference: https://github.com/scverse/anndata/blob/286bc7f207863964cb861f17c96ab24fe0cf72ac/anndata/_core/sparse_dataset.py#L230
        return arr.to_memory()
    return arr


def to_dense(arr):
    """
    Convert a sparse array to dense.

    :param arr: The array to convert.
    :type arr: np.array

    :returns: The converted array (or the original array if it was already dense).
    :rtype: np.array
    """
    if issparse(arr):
        return arr.todense()
    return arr


def to_uint8(arr, norm_along=None):
    """
    Convert an array to uint8 dtype.

    :param arr: The array to convert.
    :type arr: np.array
    :param norm_along: How to normalize the array values. By default, None. Valid values are "global", "var", "obs".
    :type norm_along: str or None

    :returns: The converted array.
    :rtype: np.array
    """
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


def sort_var_axis(adata_X, orig_var_index, full_var_index=None):
    """
    Sort the var index by performing hierarchical clustering.

    :param adata_X: The matrix to use for clustering. For example, adata.X
    :type adata_X: np.array
    :param orig_var_index: The original var index. For example, adata.var.index
    :type orig_var_index: pandas.Index
    :param full_var_index: Pass the full adata.var.index to append the var values excluded from sorting, if adata_X and orig_var_index are a subset of the full adata.X matrix. By default, None.
    :type full_var_index: pandas.Index or None

    :returns: The sorted elements of the var index.
    :rtype: list[str]
    """
    gexp_arr = to_dense(adata_X)

    # Perform hierarchical clustering along the genes axis.
    Z = scipy.cluster.hierarchy.linkage(gexp_arr.T, method="ward")

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = orig_var_index[leaf_index_list].tolist()

    if full_var_index:
        leaf_list = leaf_list + list(set(full_var_index).difference(set(leaf_list)))

    return leaf_list


def to_diamond(x, y, r):
    """
    Convert an (x, y) coordinate to a polygon (diamond) with a given radius.

    :param x: The x coordinate.
    :type x: int or float
    :param y: The y coordinate.
    :type x: int or float
    :param r: The radius.
    :type r: int or float

    :returns: The polygon vertices as an array of coordinate pairs, like [[x1, y1], [x2, y2], ...]
    :rtype: np.array
    """
    return np.array([[x, y + r], [x + r, y], [x, y - r], [x - r, y]])
