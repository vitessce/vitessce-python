import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Export data to AWS S3
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
    import boto3
    import json
    from urllib.parse import quote_plus
    from os.path import join, isfile, isdir
    from urllib.request import urlretrieve
    from anndata import read_h5ad
    import scanpy as sc

    from vitessce import (
        VitessceWidget,
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
        boto3,
        cm,
        isdir,
        isfile,
        join,
        json,
        optimize_adata,
        os,
        quote_plus,
        read_h5ad,
        urlretrieve,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Download and process data

        For this example, we need to download a dataset from the COVID-19 Cell Atlas https://www.covid19cellatlas.org/index.healthy.html#habib17.
        """
    )
    return


@app.cell
def _(isfile, join, os, read_h5ad, urlretrieve):
    adata_filepath = join("data", "habib17.processed.h5ad")
    if not isfile(adata_filepath):
        os.makedirs("data", exist_ok=True)
        urlretrieve('https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad', adata_filepath)

    adata = read_h5ad(adata_filepath)
    top_dispersion = adata.var["dispersions_norm"][
        sorted(
            range(len(adata.var["dispersions_norm"])),
            key=lambda k: adata.var["dispersions_norm"][k],
        )[-51:][0]
    ]
    adata.var["top_highly_variable"] = (
        adata.var["dispersions_norm"] > top_dispersion
    )
    return (adata,)


@app.cell
def _(VAR_CHUNK_SIZE, adata, isdir, join, optimize_adata):
    zarr_filepath = join('data', 'habib17.processed.zarr')
    if not isdir(zarr_filepath):
        adata_1 = optimize_adata(adata, obs_cols=['CellType'], obsm_keys=['X_umap'], var_cols=['top_highly_variable'], optimize_X=True)
        adata_1.write_zarr(zarr_filepath, chunks=[adata_1.shape[0], VAR_CHUNK_SIZE])
    return (zarr_filepath,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Create the Vitessce configuration
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Set up the configuration by adding the views and datasets of interest.
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, zarr_filepath):
    vc = VitessceConfig(schema_version="1.0.15", name='Habib et al', description='COVID-19 Healthy Donor Brain')
    dataset = vc.add_dataset(name='Brain').add_object(AnnDataWrapper(
            adata_path=zarr_filepath,
            obs_embedding_paths=["obsm/X_umap"],
            obs_embedding_names=["UMAP"],
            obs_set_paths=["obs/CellType"],
            obs_set_names=["Cell Type"],
            obs_feature_matrix_path="X",
            feature_filter_path="var/top_highly_variable"
    ))
    scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset)
    vc.layout((scatterplot | (cell_sets / genes)) / heatmap);
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Create a `boto3` resource with S3 credentials
        """
    )
    return


@app.cell
def _(boto3, os):
    s3 = boto3.resource(
        service_name='s3',
        aws_access_key_id=os.environ['VITESSCE_S3_ACCESS_KEY_ID'],
        aws_secret_access_key=os.environ['VITESSCE_S3_SECRET_ACCESS_KEY'],
    )
    return (s3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 5. Upload files to S3

        The `.export(to='S3')` method on the view config instance will upload all data objects to the specified bucket. Then, the processed view config will be returned as a `dict`, with the file URLs filled in, pointing to the S3 bucket files. For more information about configuring the S3 bucket so that files are accessible over the internet, visit the "Hosting Data" page of our core documentation site.
        """
    )
    return


@app.cell
def _(s3, vc):
    config_dict = vc.export(to='S3', s3=s3, bucket_name='vitessce-export-examples', prefix='test')
    return (config_dict,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 6. View on vitessce.io

        The returned view config dict can be converted to a URL, and can be used to share the interactive visualizations with colleagues.
        """
    )
    return


@app.cell
def _(config_dict, json, quote_plus):
    vitessce_url = "http://vitessce.io/?url=data:," + quote_plus(json.dumps(config_dict))
    import webbrowser
    webbrowser.open(vitessce_url)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
