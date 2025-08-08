import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of single-cell RNA seq data
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
        os,
        read_h5ad,
        urlretrieve,
    )


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
def _(isfile, join, os, urlretrieve):
    adata_filepath = join("data", "habib17.processed.h5ad")
    if not isfile(adata_filepath):
        os.makedirs("data", exist_ok=True)
        urlretrieve('https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad', adata_filepath)
    return (adata_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Load the data

        Note: this function may print a `FutureWarning`
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
        ## 3.1. Preprocess the Data For Visualization

        This dataset contains 25,587 genes.  We prepare to visualize the top 50 highly variable genes for the heatmap as ranked by dispersion norm, although one may use any boolean array filter for the heatmap.
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
        ## 3.2 Save the Data to Zarr store

        We want to convert the original `h5ad` file to a [Zarr](https://zarr.readthedocs.io/en/stable/) store, which Vitessce is able to load. We can use the `optimize_adata` function to ensure that all arrays and dataframe columns that we intend to use in our visualization are in the optimal format to be loaded by Vitessce. This function will cast arrays to numerical data types that take up less space (as long as the values allow). Note: unused arrays and columns (i.e., not specified in any of the parameters to `optimize_adata`) will not be copied into the new AnnData object.
        """
    )
    return


@app.cell
def _(VAR_CHUNK_SIZE, adata, isdir, join):
    zarr_filepath = join("data", "habib17.h5ad.zarr")
    if not isdir(zarr_filepath):
        adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
    return (zarr_filepath,)


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, zarr_filepath):
    vc = VitessceConfig(
        schema_version="1.0.17",
        name='Habib et al',
        description='COVID-19 Healthy Donor Brain'
    )

    # Add data.
    dataset = vc.add_dataset(name='Brain').add_object(AnnDataWrapper(
        adata_path=zarr_filepath,
        obs_embedding_paths=["obsm/X_umap"],
        obs_embedding_names=["UMAP"],
        obs_set_paths=["obs/CellType"],
        obs_set_names=["Cell Type"],
        obs_feature_matrix_path="X",
        initial_feature_filter_path="var/top_highly_variable",
        coordination_values={
              "obsType": 'cell',
              "featureType": 'gene',
              "featureValueType": 'expression',
        },
    )).add_object(AnnDataWrapper(
        adata_path=zarr_filepath,
        obs_feature_column_paths=["obs/percent_mito"],
        coordination_values={
            "obsType": 'cell',
            "featureType": 'qualityMetric',
            "featureValueType": 'value',
        }
    ))

    # Add views.
    scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    scatterplot_2 = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
    histogram = vc.add_view(cm.FEATURE_VALUE_HISTOGRAM, dataset=dataset)

    # Link views.

    # Color one of the two scatterplots by the percent_mito quality metric.
    # Also use this quality metric for the histogram values.
    vc.link_views_by_dict([histogram, scatterplot_2], {
        "obsType": 'cell',
        "featureType": 'qualityMetric',
        "featureValueType": 'value',
        "featureSelection": ["percent_mito"],
        "obsColorEncoding": "geneSelection",
    }, meta=False)

    # Synchronize the zooming and panning of the two scatterplots
    vc.link_views_by_dict([scatterplot, scatterplot_2], {
        "embeddingZoom": None,
        "embeddingTargetX": None,
        "embeddingTargetY": None,
    }, meta=False)

    # Define the layout.
    vc.layout((scatterplot | (cell_sets / genes)) / (scatterplot_2 | histogram));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 5. Create the widget

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
