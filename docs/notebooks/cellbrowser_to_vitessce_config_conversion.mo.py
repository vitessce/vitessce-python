import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Load UCSC Cell Browser project in Vitessce
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        This notebook shows you how to use the `convert_cell_browser_project_to_anndata` function, which allows you to take an existing project, published in https://cells.ucsc.edu/ and:
        1. Convert it into the AnnData format that is supported by Vitessce
        2. Save the AnnData object as a Zarr store
        3. Configure Vitessce with the AnnData-Zarr store
        4. Render a Vitessce widget based on the config (step 3) directly in the notebook.

        The dataset that you choose to convert needs to be a valid UCSC Cell Browser "project", accessible from https://cells.ucsc.edu/, with a configuration available in https://github.com/ucscGenomeBrowser/cellbrowser-confs

        The `convert_cell_browser_project_to_anndata` function takes the name of that project as an input. For example, to convert this project, https://cells.ucsc.edu/?ds=adultPancreas, you will neeed to pass `"adultPancreas"` as the project name.
        """
    )
    return


@app.cell
def _():
    import os
    import json
    from os.path import join
    from vitessce import (
        convert_cell_browser_project_to_anndata,
        AnnDataWrapper,
        VitessceConfig,
    )
    from vitessce.data_utils import VAR_CHUNK_SIZE
    return (
        AnnDataWrapper,
        VAR_CHUNK_SIZE,
        VitessceConfig,
        convert_cell_browser_project_to_anndata,
        join,
        os,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Convert UCSC Cell Browser project to a format that is supported by Vitessce
        #### Output:
        An AnnData object

        """
    )
    return


@app.cell
def _():
    ## 3. Convert UCSC Cell Browser project to a Vitessce view config
    return


@app.cell
def _(convert_cell_browser_project_to_anndata):
    # Example run, coverting "adultPancreas" project:
    adata = convert_cell_browser_project_to_anndata(project_name="adultPancreas", keep_only_marker_genes=True)
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Save the AnnData object as a Zarr store
        """
    )
    return


@app.cell
def _(VAR_CHUNK_SIZE, adata, join, os):
    zarr_filepath = join("data", "out.adata.zarr")
    os.makedirs(os.path.dirname(zarr_filepath), exist_ok=True)
    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
    return (zarr_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Configure Vitessce with the AnnData-Zarr store
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, zarr_filepath):
    anndata_wrapper_inst = AnnDataWrapper(
        adata_path=zarr_filepath,
        obs_feature_matrix_path="X",
        obs_embedding_paths=["obsm/X_tsne"],
        obs_embedding_names=["t-SNE"],
        obs_set_paths=["obs/cluster", "obs/age"],
        obs_set_names=["cluster", "age"],
    )
    vc = VitessceConfig(schema_version="1.0.15", name="Vitessce configuration for CellBrowser project adultPancreas")
    anndata_wrapper_inst.auto_view_config(vc)
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Render the Vitessce widget
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
