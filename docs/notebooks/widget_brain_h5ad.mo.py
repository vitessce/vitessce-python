import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of single-cell RNA seq data from H5AD file
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
    import json

    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
    )
    from vitessce.data_utils import (
        generate_h5ad_ref_spec
    )
    return (
        AnnDataWrapper,
        VitessceConfig,
        cm,
        generate_h5ad_ref_spec,
        isfile,
        join,
        json,
        os,
        urlretrieve,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 0. Download data
        """
    )
    return


@app.cell
def _():
    h5_url = "https://datasets.cellxgene.cziscience.com/84df8fa1-ab53-43c9-a439-95dcb9148265.h5ad"
    return (h5_url,)


@app.cell
def _(h5_url, isfile, join, os, urlretrieve):
    adata_filepath = join("data", "84df8fa1-ab53-43c9-a439-95dcb9148265.h5ad")
    if not isfile(adata_filepath):
        os.makedirs("data", exist_ok=True)
        urlretrieve(h5_url, adata_filepath)
    return (adata_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Create a Reference Spec JSON file for the H5AD file

        In order for Vitessce to load H5AD files, we also need to provide a corresponding [Reference Spec](https://fsspec.github.io/kerchunk/spec.html) JSON file which contains mappings between AnnData object keys and the byte offsets at which those AnnData object values begin within the H5AD file binary contents.
        """
    )
    return


@app.cell
def _(generate_h5ad_ref_spec, h5_url, isfile, join, json):
    json_filepath = join("data", "84df8fa1-ab53-43c9-a439-95dcb9148265.h5ad.reference.json")
    if not isfile(json_filepath):
        ref_dict = generate_h5ad_ref_spec(h5_url)
        with open(json_filepath, "w") as f:
            json.dump(ref_dict, f)
    return (json_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Create the Vitessce widget configuration

        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, adata_filepath, cm, json_filepath):
    vc = VitessceConfig(schema_version="1.0.17", name='Nakshatri et al', description='snRNA-seq analyses of breast tissues of healthy women of diverse genetic ancestry')

    dataset = vc.add_dataset(name='84df8fa1').add_object(AnnDataWrapper(
            adata_path=adata_filepath,
            ref_path=json_filepath, # We specify paths to both the H5AD and JSON files
            obs_embedding_paths=["obsm/X_wnn.umap"],
            obs_embedding_names=["UMAP"],
            obs_set_paths=["obs/cell_type"],
            obs_set_names=["Cell Type"],
            obs_feature_matrix_path="X",
        )
    )

    scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    cell_set_sizes = vc.add_view(cm.OBS_SET_SIZES, dataset=dataset)
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)

    vc.layout((scatterplot | cell_sets) / (cell_set_sizes | genes));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Create the widget
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
