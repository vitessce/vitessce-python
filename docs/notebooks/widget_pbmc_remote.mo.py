import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of 3k PBMC reference from Remote Zarr Store
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
    from os.path import join
    from urllib.request import urlretrieve
    from anndata import read_h5ad
    import scanpy as sc

    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
    )
    return AnnDataWrapper, VitessceConfig, cm


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Set the URL for the Remote Dataset

        For this example, we already have uploaded the `pbmc3k` dataset as a zarr store from the [scanpy docs](https://scanpy.readthedocs.io/en/stable/api/scanpy.datasets.pbmc3k.html) to the cloud.
        """
    )
    return


@app.cell
def _():
    url = 'https://storage.googleapis.com/vitessce-demo-data/anndata-test/pbmc3k_processed.zarr/'
    return (url,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Create a Vitessce view config

        Define the data and views you would like to include in the widget.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, url):
    vc = VitessceConfig(schema_version="1.0.15", name='PBMC Reference')
    dataset = vc.add_dataset(name='PBMC 3k').add_object(AnnDataWrapper(adata_url=url, obs_set_paths=["obs/louvain"], obs_set_names=["Louvain"], obs_embedding_paths=["obsm/X_umap", "obsm/X_pca"], obs_embedding_names=["UMAP", "PCA"], obs_feature_matrix_path="X"))

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
        ## 4. Create the Vitessce widget
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        A widget can be created with the `.widget()` method on the config instance. Here, the `proxy=True` parameter allows this widget to be used in a cloud notebook environment, such as Binder.
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
