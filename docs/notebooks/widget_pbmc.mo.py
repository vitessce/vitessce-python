import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Vitessce Widget Tutorial
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of 3k PBMC reference
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Import dependencies

        We need to import the classes and functions that we will be using from the corresponding packages.
        """
    )
    return


@app.cell
def _():
    import os
    from os.path import join, isfile, isdir
    from urllib.request import urlretrieve
    from anndata import read_h5ad
    import scanpy as sc

    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
    )
    from vitessce.data_utils import (
        optimize_adata,
        VAR_CHUNK_SIZE,
    )
    return (
        AnnDataWrapper,
        VAR_CHUNK_SIZE,
        VitessceConfig,
        cm,
        isdir,
        isfile,
        join,
        optimize_adata,
        os,
        read_h5ad,
        urlretrieve,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Download the dataset

        Download `pbmc3k_final.h5ad` from https://seurat.nygenome.org/pbmc3k_final.h5ad
        """
    )
    return


@app.cell
def _(isfile, join, os, urlretrieve):
    adata_filepath = join("data", "pbmc3k_final.h5ad")
    if not isfile(adata_filepath):
        os.makedirs("data", exist_ok=True)
        urlretrieve('https://seurat.nygenome.org/pbmc3k_final.h5ad', adata_filepath)
    return (adata_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Load the dataset

        Load the dataset using AnnData's `read_h5ad` function.
        """
    )
    return


@app.cell
def _(adata_filepath, read_h5ad):
    adata = read_h5ad(adata_filepath)
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3.1 Save the AnnData object to Zarr
        """
    )
    return


@app.cell
def _(VAR_CHUNK_SIZE, adata, isdir, join, optimize_adata):
    zarr_filepath = join('data', 'pbmc3k_final.zarr')
    if not isdir(zarr_filepath):
        adata_1 = optimize_adata(adata, obs_cols=['leiden'], obsm_keys=['X_umap', 'X_pca'], optimize_X=True)
        adata_1.write_zarr(zarr_filepath, chunks=[adata_1.shape[0], VAR_CHUNK_SIZE])
    return (zarr_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Create a Vitessce view config

        Define the data and views you would like to include in the widget.

        For more details about how to configure data depending on where the files are located relative to the notebook execution, see https://python-docs.vitessce.io/data_options.html.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, zarr_filepath):
    vc = VitessceConfig(schema_version="1.0.15", name='PBMC Reference')
    dataset = vc.add_dataset(name='PBMC 3k').add_object(AnnDataWrapper(
        adata_store=zarr_filepath,
        obs_set_paths=["obs/leiden"],
        obs_set_names=["Leiden"],
        obs_embedding_paths=["obsm/X_umap", "obsm/X_pca"],
        obs_embedding_names=["UMAP", "PCA"],
        obs_feature_matrix_path="X"
    ))

    umap = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    pca = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="PCA")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset)

    vc.layout((umap / pca) | ((cell_sets | genes) / heatmap));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 5. Create the Vitessce widget
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
