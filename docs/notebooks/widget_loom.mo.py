import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of a Loom file
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Import dependencies
        """
    )
    return


@app.cell
def _():
    import os
    from os.path import join, isfile, isdir
    from urllib.request import urlretrieve
    from anndata import read_loom
    import numpy as np

    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
    )
    from vitessce.data_utils import (
        optimize_adata,
        to_diamond,
        VAR_CHUNK_SIZE,
    )
    return (
        AnnDataWrapper,
        VAR_CHUNK_SIZE,
        VitessceConfig,
        cm,
        ct,
        isdir,
        isfile,
        join,
        np,
        optimize_adata,
        os,
        read_loom,
        to_diamond,
        urlretrieve,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Download data

        Download `osmFISH_SScortex_mouse_all_cells.loom` from http://loom.linnarssonlab.org/.
        """
    )
    return


@app.cell
def _(isfile, join, os, urlretrieve):
    loom_filepath = join("data", "osmFISH_SScortex_mouse_all_cells.loom")
    if not isfile(loom_filepath):
        os.makedirs("data", exist_ok=True)
        urlretrieve('http://loom.linnarssonlab.org/clone/osmFISH/osmFISH_SScortex_mouse_all_cells.loom', loom_filepath)
    return (loom_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Open Loom file with AnnData's read_loom
        """
    )
    return


@app.cell
def _(loom_filepath, read_loom):
    adata = read_loom(loom_filepath, obsm_names={"tSNE": ["_tSNE_1", "_tSNE_2"], "spatial": ["X", "Y"]})
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Generate pseudo-segmentations as diamond-shaped polygons centered on the spatial coordinate of each cell, and store in `adata.obsm["segmentations"]`
        """
    )
    return


@app.cell
def _(adata, np, to_diamond):
    num_cells = adata.obs.shape[0]
    adata.obsm["segmentations"] = np.zeros((num_cells, 4, 2))
    radius = 100
    for i in range(num_cells):
        adata.obsm["segmentations"][i, :, :] = to_diamond(adata.obsm['spatial'][i, 0], adata.obsm['spatial'][i, 1], radius)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Save the AnnData object to a Zarr store:
        """
    )
    return


@app.cell
def _(VAR_CHUNK_SIZE, adata, isdir, join, optimize_adata):
    zarr_filepath = join('data', 'osmFISH_SScortex_mouse_all_cells.zarr')
    if not isdir(zarr_filepath) or True:
        adata_1 = optimize_adata(adata, obs_cols=['ClusterName'], obsm_keys=['tSNE', 'spatial', 'segmentations'], optimize_X=True)
        adata_1.write_zarr(zarr_filepath, chunks=[adata_1.shape[0], VAR_CHUNK_SIZE])
    return (zarr_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Configure Vitessce

        Create a Vitessce view config.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, ct, zarr_filepath):
    vc = VitessceConfig(schema_version="1.0.15", name='Loom Example', description='osmFISH dataset of the mouse cortex including all cells')
    w = AnnDataWrapper(adata_path=zarr_filepath, obs_set_paths=["obs/ClusterName"], obs_set_names=["Clusters"], obs_locations_path="obsm/spatial", obs_segmentations_path="obsm/segmentations", obs_embedding_paths=["obsm/tSNE"])
    dataset = vc.add_dataset(name='SScortex').add_object(w)

    tsne = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="tSNE")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    spatial = vc.add_view(cm.SPATIAL, dataset=dataset)

    spatial_segmentation_layer_value = {
      "opacity": 1,
      "radius": 0,
      "visible": True,
      "stroked": False
    }

    vc.link_views([spatial], [ct.SPATIAL_ZOOM, ct.SPATIAL_TARGET_X, ct.SPATIAL_TARGET_Y, ct.SPATIAL_SEGMENTATION_LAYER], [-6.43, 10417.69, 24885.55, spatial_segmentation_layer_value])
    vc.layout(spatial | (tsne / cell_sets));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 5. Render the widget
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        A widget can be created with the `.widget()` method on the config instance.
        """
    )
    return


@app.cell
def _(vc):
    vw = vc.widget()
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
