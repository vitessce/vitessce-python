import argparse
from anndata import read_loom
import pickle
import numpy as np
import scanpy as sc
import geopandas as gpd
from shapely.geometry import Polygon
from vitessce.data_utils import (
    to_uint8,
    sort_var_axis,
    optimize_adata,
)

# Taken from http://linnarssonlab.org/osmFISH/clusters/
LOOKUP = {
    "Astrocyte Gfap": "Astrocyte",
    "Astrocyte Mfge8": "Astrocyte",
    "C. Plexus": "Ventricle",
    "Endothelial 1": "Vasculature",
    "Endothelial": "Vasculature",
    "Ependymal": "Ventricle",
    "Hippocampus": "Excitatory neurons",
    "Inhibitory CP": "Inhibitory neurons",
    "Inhibitory Cnr1": "Inhibitory neurons",
    "Inhibitory Crhbp": "Inhibitory neurons",
    "Inhibitory IC": "Inhibitory neurons",
    "Inhibitory Kcnip2": "Inhibitory neurons",
    "Inhibitory Pthlh": "Inhibitory neurons",
    "Inhibitory Vip": "Inhibitory neurons",
    "Microglia": "Brain immune",
    "Oligodendrocyte COP": "Oligodendrocytes",
    "Oligodendrocyte MF": "Oligodendrocytes",
    "Oligodendrocyte Mature": "Oligodendrocytes",
    "Oligodendrocyte NF": "Oligodendrocytes",
    "Oligodendrocyte Precursor cells": "Oligodendrocytes",
    "Pericytes": "Vasculature",
    "Perivascular Macrophages": "Brain immune",
    "Pyramidal Cpne5": "Excitatory neurons",
    "Pyramidal Kcnip2": "Excitatory neurons",
    "Pyramidal L2-3 L5": "Excitatory neurons",
    "Pyramidal L2-3": "Excitatory neurons",
    "Pyramidal L3-4": "Excitatory neurons",
    "Pyramidal L5": "Excitatory neurons",
    "Pyramidal L6": "Excitatory neurons",
    "Vascular Smooth Muscle": "Vasculature",
    "pyramidal L4": "Excitatory neurons",
}

# Copied from https://github.com/vitessce/vitessce-data/blob/6726dae/python/cell_reader.py#L17


def octagon(poly):
    """
    Returns a bounding octagon.
    >>> square = np.array([[0,0], [0,1], [1,1], [1,0]])
    >>> octagon(square)
    [[0, 0], [0, 1], [0, 1], [1, 1], [1, 1], [1, 0], [1, 0], [0, 0]]
    >>> triangle = np.array([[1,0], [0,2], [2,3]])
    >>> octagon(triangle)
    [[0, 1], [0, 2], [1, 3], [2, 3], [2, 3], [2, 1], [1, 0], [1, 0]]
    >>> type(octagon(triangle)[0][0])
    <class 'int'>
    """
    # SciPy has ConvexHull, but it segfaulted on me: perhaps
    #   https://github.com/scipy/scipy/issues/9751
    # ... and even if I fixed it locally,
    # not sure we want that fragility.
    #
    # Also: The goal is really just to get a simplified shape...
    # A convex hull is too precise in some ways,
    # while in others it falls short, ie, concavities.
    #
    # I kind-of like the obvious artificiality of an octagon.

    # Was unsigned, and substraction causes underflow.
    poly_as_int = poly
    min_x = np.min(poly_as_int[:, [0]])
    max_x = np.max(poly_as_int[:, [0]])
    min_y = np.min(poly_as_int[:, [1]])
    max_y = np.max(poly_as_int[:, [1]])

    summed = np.sum(poly_as_int, axis=1)
    diffed = np.diff(poly_as_int, axis=1)

    min_sum = np.min(summed)
    max_sum = np.max(summed)
    min_diff = np.min(diffed)
    max_diff = np.max(diffed)

    return Polygon(
        [
            [min_x, min_sum - min_x],
            [min_x, max_diff + min_x],  # ^ Left
            [max_y - max_diff, max_y],
            [max_sum - max_y, max_y],  # ^ Botton
            [max_x, max_sum - max_x],
            [max_x, min_diff + max_x],  # ^ Right
            [min_y - min_diff, min_y],
            [min_sum - min_y, min_y],  # ^ Top
        ]
    )


def convert_to_cells_h5ad_zarr(args):
    adata = read_loom(
        args.input_loom,
        obsm_mapping={
            "X_tsne": ["_tSNE_1", "_tSNE_2"],
            "X_spatial": ["X", "Y"],
        },
    )

    adata.obs["is_excluded"] = adata.obs["ClusterName"] == "Excluded"
    adata = adata[~adata.obs["is_excluded"], :].copy()

    # Reorder the genes axis after hierarchical clustering.
    leaf_list = sort_var_axis(adata.X, adata.var.index.values)
    adata = adata[:, leaf_list].copy()

    # Store expression matrix as uint8
    adata.layers["X_uint8"] = to_uint8(adata.X, norm_along="var")

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    adata.obs = adata.obs.rename(columns={"ClusterName": "Subcluster"})
    adata.obs["Cluster"] = adata.obs["Subcluster"].apply(lambda val: LOOKUP[val])

    adata.obsm["X_centroid"] = np.zeros((adata.shape[0], 2))
    adata.obsm["X_segmentations"] = np.zeros((adata.shape[0], 8, 2))

    with open(args.input_pkl, "rb") as f_in:
        segmentation = pickle.load(f_in)

        cell_ids = adata.obs.index.values.tolist()

        geoseries = gpd.GeoSeries(
            index=cell_ids,
            data=[Polygon(segmentation[cell_id]) for cell_id in cell_ids],
        )
        geoseries_simple = geoseries.apply(
            lambda p: octagon(np.array(p.exterior.coords))
        )
        centroids = geoseries_simple.centroid

        adata.obsm["X_centroid"][:, 0] = centroids.x
        adata.obsm["X_centroid"][:, 1] = centroids.y

        for cell_id in cell_ids:
            i = cell_ids.index(cell_id)
            adata.obsm["X_segmentations"][i, :, :] = np.array(
                geoseries_simple.loc[cell_id].exterior.coords
            )[:8, :]

    adata = optimize_adata(
        adata,
        obs_cols=["Cluster", "Subcluster", "Region"],
        var_cols=["Fluorophore", "Hybridization"],
        obsm_keys=[
            "X_tsne",
            "X_umap",
            "X_pca",
            "X_spatial",
            "X_centroid",
            "X_segmentations",
        ],
        layer_keys=["X_uint8"],
    )

    adata.write_zarr(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-il", "--input_loom", type=str, required=True, help="Input Loom file"
    )
    parser.add_argument(
        "-ip", "--input_pkl", type=str, required=True, help="Input pickle file"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output AnnData Zarr store",
    )
    args = parser.parse_args()
    convert_to_cells_h5ad_zarr(args)
