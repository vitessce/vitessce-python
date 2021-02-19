# This script has been adapted from https://github.com/vitessce/vitessce-data/blob/6726dae/python/cell_reader.py
import argparse

import json
import pickle
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc

from anndata import read_loom

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

def octagon(poly):
    '''
    Returns a bounding octagon.
    >>> square = np.array([[0,0], [0,1], [1,1], [1,0]])
    >>> octagon(square)
    [[0, 0], [0, 1], [0, 1], [1, 1], [1, 1], [1, 0], [1, 0], [0, 0]]
    >>> triangle = np.array([[1,0], [0,2], [2,3]])
    >>> octagon(triangle)
    [[0, 1], [0, 2], [1, 3], [2, 3], [2, 3], [2, 1], [1, 0], [1, 0]]
    >>> type(octagon(triangle)[0][0])
    <class 'int'>
    '''
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
    poly_as_int = poly.astype('int')
    min_x = int(np.min(poly_as_int[:, [0]]))
    max_x = int(np.max(poly_as_int[:, [0]]))
    min_y = int(np.min(poly_as_int[:, [1]]))
    max_y = int(np.max(poly_as_int[:, [1]]))

    summed = np.sum(poly_as_int, axis=1)
    diffed = np.diff(poly_as_int, axis=1)

    min_sum = int(np.min(summed))
    max_sum = int(np.max(summed))
    min_diff = int(np.min(diffed))
    max_diff = int(np.max(diffed))

    return [
        [min_x, min_sum - min_x],
        [min_x, max_diff + min_x],  # ^ Left
        [max_y - max_diff, max_y],
        [max_sum - max_y, max_y],  # ^ Botton
        [max_x, max_sum - max_x],
        [max_x, min_diff + max_x],  # ^ Right
        [min_y - min_diff, min_y],
        [min_sum - min_y, min_y]  # ^ Top
    ]


def mean_coord(coords):
    '''
    The xy values in the Linnarsson data are not good:
    They take a different corner as the origin.
    So... we find the center of our polygon instead.
    >>> mean_coord([[1,2], [3,4], [5,6]])
    [3, 4]
    '''
    return [int(x) for x in np.mean(coords, axis=0).tolist()]

def get_cells_dict(adata, segmentation, round_xy=True):
    cells = {}

    for i, cell_id in enumerate(adata.obs.index.values):
        subcluster = adata.obs.at[cell_id, 'ClusterName']
        tsne_coords = adata.obsm['X_tsne'][i, :].astype(float)
        pca_coords = adata.obsm['X_pca'][i, :].astype(float)
        umap_coords = adata.obsm['X_umap'][i, :].astype(float)

        cells[cell_id] = {
            'mappings': {
                't-SNE': [tsne_coords[0], tsne_coords[1]],
                'PCA': [pca_coords[0], pca_coords[1]],
                'UMAP': [umap_coords[0], umap_coords[1]],
            },
            'factors': {
                'subcluster': subcluster,
                'cluster': LOOKUP[subcluster]
            }
        }
    
    for cell_id, poly in segmentation.items():
        if cell_id in cells:
            simple_poly = octagon(poly)
            xy = mean_coord(simple_poly)
            if round_xy:
                xy = [ int(z) for z in xy ]
            cells[cell_id]['poly'] = simple_poly
            cells[cell_id]['xy'] = xy
    
    return cells

def get_cell_sets_dict(adata):
    return {}

def get_exp_matrix_dict(adata):
    return {}

def get_genes_dict(adata):
    return {}

def get_factors_dict(adata):
    return {}

def get_neighborhoods_dict(adata):
    return {}

def convert_loom_to_cells_json(args):
    # Load raw files.
    adata = read_loom(args.input_loom, obsm_names={"X_tsne": ["_tSNE_1", "_tSNE_2"], "spatial": ["X", "Y"]})
    with open(args.input_pkl, 'rb') as f:
        segmentation = pickle.load(f)

    # Pre-process data.
    adata = adata[adata.obs['ClusterName'] != 'Excluded', :]
    sc.pp.neighbors(adata)
    sc.tl.pca(adata, random_state=2445)
    sc.tl.umap(adata, random_state=2445)


    # Get dict representations and save to JSON output files.
    cells_dict = get_cells_dict(adata, segmentation)
    json.dump(cells_dict, args.output_cells)

    cell_sets_dict = get_cell_sets_dict(adata)
    json.dump(cell_sets_dict, args.output_cell_sets)

    exp_matrix_dict = get_exp_matrix_dict(adata)
    json.dump(exp_matrix_dict, args.output_exp_matrix)

    genes_dict = get_genes_dict(adata)
    json.dump(genes_dict, args.output_genes)

    factors_dict = get_factors_dict(adata)
    json.dump(factors_dict, args.output_factors)
    
    neighborhoods_dict = get_neighborhoods_dict(adata)
    json.dump(neighborhoods_dict, args.output_neighborhoods)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-il', '--input_loom', type=str, required=True, help='Input Loom file')
    parser.add_argument('-ip', '--input_pkl', type=str, required=True, help='Input pickle file containing cell segmentations')
    
    parser.add_argument('-oc', '--output_cells', type=argparse.FileType('w'), required=True, help='Output cells.json file')
    parser.add_argument('-ocs', '--output_cell_sets', type=argparse.FileType('w'), required=True, help='Output cell-sets.json file')
    parser.add_argument('-oem', '--output_exp_matrix', type=argparse.FileType('w'), required=True, help='Output clusters.json file')
    parser.add_argument('-og', '--output_genes', type=argparse.FileType('w'), required=True, help='Output genes.json file')
    parser.add_argument('-on', '--output_neighborhoods', type=argparse.FileType('w'), required=True, help='Output neighborhoods.json file')
    parser.add_argument('-of', '--output_factors', type=argparse.FileType('w'), required=True, help='Output factors.json file')

    args = parser.parse_args()
    convert_loom_to_cells_json(args)
