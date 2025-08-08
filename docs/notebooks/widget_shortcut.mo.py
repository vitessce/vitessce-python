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
        # The from_object shortcut
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Import dependencies

        Import the functions and classes that we will be using.
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
    return AnnDataWrapper, VitessceConfig, join, os, read_h5ad, urlretrieve


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Download the data

        For this example, we need to download a dataset from the COVID-19 Cell Atlas https://www.covid19cellatlas.org/index.healthy.html#habib17.
        """
    )
    return


@app.cell
def _(join, os, urlretrieve):
    os.makedirs("data", exist_ok=True)
    adata_filepath = join("data", "habib17.processed.h5ad")
    urlretrieve('https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad', adata_filepath)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Load the data
        """
    )
    return


@app.cell
def _(join, read_h5ad):
    adata = read_h5ad(join("data", "habib17.processed.h5ad"))
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3.1. Preprocess the Data For Visualization
        """
    )
    return


@app.cell
def _(adata):
    top_dispersion = adata.var["dispersions_norm"][
        sorted(
            range(len(adata.var["dispersions_norm"])),
            key=lambda k: adata.var["dispersions_norm"][k],
        )[-51:][0]
    ]
    adata.var["top_highly_variable"] = (
        adata.var["dispersions_norm"] > top_dispersion
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        With one line of code, you may create a Vitessce widget based on an automatically inferred configuration.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, adata):
    vw = VitessceConfig.from_object(AnnDataWrapper(
            adata,
            obs_embedding_paths=["obsm/X_umap"],
            obs_embedding_names=["UMAP"],
            obs_set_paths=["obs/CellType"],
            obs_set_names=["Cell Type"],
            obs_feature_matrix_path="X",
            feature_filter_path="var/top_highly_variable"
    ), schema_version="1.0.15").widget(height=800)
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
