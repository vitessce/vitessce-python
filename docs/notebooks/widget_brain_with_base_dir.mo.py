import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Configure relative to a base_dir
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

    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
        BASE_URL_PLACEHOLDER,
    )
    from vitessce.data_utils import (
        optimize_adata,
        VAR_CHUNK_SIZE,
    )
    return (
        AnnDataWrapper,
        BASE_URL_PLACEHOLDER,
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
        ## 2. Define a `base_dir`

        We will define a `base_dir` inside which our data will live. We will provide this to `VitessceConfig` in order to construct a configuration that contains URL paths relative to this directory.
        """
    )
    return


@app.cell
def _():
    BASE_DIR = "data"
    return (BASE_DIR,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Download the data

        For this example, we need to download a dataset from the COVID-19 Cell Atlas https://www.covid19cellatlas.org/index.healthy.html#habib17.
        """
    )
    return


@app.cell
def _(BASE_DIR, isfile, join, os, urlretrieve):
    adata_relative_filepath = "habib17.processed.h5ad" # Relative to BASE_DIR
    adata_filepath = join(BASE_DIR, adata_relative_filepath)
    if not isfile(adata_filepath):
        os.makedirs(BASE_DIR, exist_ok=True)
        urlretrieve('https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad', adata_filepath)
    return (adata_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Load the data

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
        ## 4.1. Preprocess the Data For Visualization

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
        ## 4.2 Save the Data to Zarr store

        We want to convert the original `h5ad` file to a [Zarr](https://zarr.readthedocs.io/en/stable/) store, which Vitessce is able to load. We can use the `optimize_adata` function to ensure that all arrays and dataframe columns that we intend to use in our visualization are in the optimal format to be loaded by Vitessce. This function will cast arrays to numerical data types that take up less space (as long as the values allow). Note: unused arrays and columns (i.e., not specified in any of the parameters to `optimize_adata`) will not be copied into the new AnnData object.
        """
    )
    return


@app.cell
def _(BASE_DIR, VAR_CHUNK_SIZE, adata, isdir, join, optimize_adata):
    zarr_relative_filepath = 'habib17.processed.zarr'
    zarr_filepath = join(BASE_DIR, zarr_relative_filepath)
    if not isdir(zarr_filepath):
        adata_1 = optimize_adata(adata, obs_cols=['CellType'], obsm_keys=['X_umap'], optimize_X=True, var_cols=['top_highly_variable'])
        adata_1.write_zarr(zarr_filepath, chunks=[adata_1.shape[0], VAR_CHUNK_SIZE])
    return (zarr_relative_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 5. Create the Vitessce widget configuration

        Vitessce needs to know which pieces of data we are interested in visualizing, the visualization types we would like to use, and how we want to coordinate (or link) the views.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 5.1. Instantiate a `VitessceConfig` object

        Use the `VitessceConfig` constructor to create an instance. In this case, we want to construct our configuration using local data that is relative to a particular directory, so we provide the `base_dir` parameter.

        Note: This `base_dir` parameter is optional. When it is omitted, local data paths are assumed to be relative to the current working directory.
        """
    )
    return


@app.cell
def _(BASE_DIR, VitessceConfig):
    vc = VitessceConfig(schema_version="1.0.15", name='Habib et al', base_dir=BASE_DIR)
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 5.2. Add a dataset to the `VitessceConfig` instance

        In Vitessce, a dataset is a container for one file per data type. The `.add_dataset(name)` method on the `vc` instance sets up and returns a new dataset instance.

        Then, we can call the dataset's `.add_object(wrapper_object)` method to attach a "data wrapper" instance to our new dataset. For example, the `AnnDataWrapper` helps to configure AnnData Zarr stores for use in the Vitessce configuration.

        Dataset wrapper classes may require additional parameters to resolve ambiguities. For instance, `AnnData` objects may store multiple clusterings or cell type annotation columns in the `adata.obs` dataframe. We can use the parameter `obs_set_paths` to tell Vitessce that certain columns of the `obs` dataframe correspond to cell type annotations or cell clusterings.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, vc, zarr_relative_filepath):
    dataset = vc.add_dataset(name='Brain').add_object(AnnDataWrapper(
            adata_path=zarr_relative_filepath, # Relative to BASE_DIR (because we specified base_dir in the VitessceConfig constructor)
            obs_embedding_paths=["obsm/X_umap"],
            obs_embedding_names=["UMAP"],
            obs_set_paths=["obs/CellType"],
            obs_set_names=["Cell Type"],
            obs_feature_matrix_path="X",
            initial_feature_filter_path="var/top_highly_variable"
        )
    )
    return (dataset,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 5.3. Add visualizations to the `VitessceConfig` instance

        Now that we have added a dataset, we can configure visualizations. The `.add_view` method adds a view (i.e. visualization or controller component) to the configuration.

        The `Component` enum class (which we have imported as `cm` here) can be used to fill in the `component_type` parameter.

        For convenience, the `SCATTERPLOT` component type takes the extra `mapping` keyword argument, which specifies which embedding should be used for mapping cells to (x,y) points on the plot.
        """
    )
    return


@app.cell
def _(cm, dataset, vc):
    scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset)
    return cell_sets, genes, heatmap, scatterplot


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### 5.4. Define the visualization layout

        The `vc.layout(view_concat)` method allows us to specify how our views will be arranged in the layout grid in the widget. The `|` and `/` characters are magic syntax for `hconcat(v1, v2)` and `vconcat(v1, v2)`, respectively.
        """
    )
    return


@app.cell
def _(cell_sets, genes, heatmap, scatterplot, vc):
    vc.layout((scatterplot | cell_sets) / (heatmap | genes));
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 6. Create the widget

        The `vc.widget()` method returns the configured widget instance.
        """
    )
    return


@app.cell
def _(vc):
    vw = vc.widget()
    vw
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 7. Check the URLs in the configuration

        We can check that the data URLs in the configuration respected the specified `base_dir`.
        """
    )
    return


@app.cell
def _(BASE_URL_PLACEHOLDER, vc):
    config_dict = vc.to_dict(base_url=BASE_URL_PLACEHOLDER)
    config_dict
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
