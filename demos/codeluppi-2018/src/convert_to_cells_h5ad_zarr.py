import argparse
from anndata import AnnData, read_loom
import pickle
import json
import numpy as np
import pandas as pd
import scanpy as sc
import geopandas as gpd
import dask_geopandas as dgpd
from dask.diagnostics import ProgressBar
from tqdm import tqdm
from shapely.geometry import Point, Polygon, shape
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
  "pyramidal L4": "Excitatory neurons"
}

def convert_to_cells_h5ad_zarr(args):
    adata = read_loom(
        args.input_loom,
        obsm_mapping={
            "X_tsne": ["_tSNE_1", "_tSNE_2"],
            "X_spatial": ["X", "Y"],
        }
    )

    adata.obs['is_excluded'] = adata.obs['ClusterName'] == 'Excluded'
    adata = adata[~adata.obs['is_excluded'], :].copy()

    # Reorder the genes axis after hierarchical clustering.
    leaf_list = sort_var_axis(adata.X, adata.var.index.values)
    adata = adata[:, leaf_list].copy()

    # Store expression matrix as uint8
    adata.layers['X_uint8'] = to_uint8(adata.X, norm_along="var")

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    adata.obs['cell_type'] = adata.obs['ClusterName'].apply(lambda val: LOOKUP[val])
    adata.obsm['X_centroid'] = np.zeros((adata.shape[0], 2))

    segmentation_simple = {}

    with open(args.input_pkl, "rb") as f_in:
        segmentation = pickle.load(f_in)

        cell_ids = adata.obs.index.values.tolist()

        print("geoseries")
        geoseries = gpd.GeoSeries([ Polygon(segmentation[cell_id]) for cell_id in cell_ids ])
        print("dgs")
        dgs = dgpd.from_geopandas(geoseries, npartitions=500)
        print("geoseries_simple")
        with ProgressBar():
            geoseries_simple = dgs.simplify(500.0, preserve_topology=True).compute()
        print("centroids")
        centroids = geoseries_simple.centroid

        adata.obsm['X_centroid'][:, 0] = centroids.x
        adata.obsm['X_centroid'][:, 1] = centroids.y

        for cell_id in cell_ids:
            segmentation_simple[cell_id] = geoseries_simple[cell_id].exterior.coords
           
    with open(args.output_segmentations, "w") as f_out:
        json.dump(segmentation_simple, f_out)
    
    adata = optimize_adata(
        adata,
        obs_cols=["ClusterName", "cell_type", "Region"],
        var_cols=["Flourophore", "Hybridization"],
        obsm_keys=["X_tsne", "X_umap", "X_pca", "X_spatial", "X_centroid"],
        layer_keys=["X_uint8"],
    )
    
    adata.write_zarr(args.output_cells)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-il',
        '--input_loom',
        type=str,
        required=True,
        help='Input Loom file'
    )
    parser.add_argument(
        '-ip',
        '--input_pkl',
        type=str,
        required=True,
        help='Input pickle file'
    )
    parser.add_argument(
        '-oc',
        '--output_cells',
        type=str,
        required=True,
        help='Output AnnData Zarr store'
    )
    parser.add_argument(
        '-os',
        '--output_segmentations',
        type=str,
        required=True,
        help='Output JSON file'
    )
    args = parser.parse_args()
    convert_to_cells_h5ad_zarr(args)
