import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of genomic profiles
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
    from vitessce import (
        VitessceConfig,
        ViewType as vt,
        CoordinationType as ct,
        AnnDataWrapper,
        MultivecZarrWrapper,
    )
    from vitessce.data_utils import (
        adata_to_multivec_zarr,
    )
    from os.path import join
    from scipy.io import mmread
    import pandas as pd
    import numpy as np
    from anndata import AnnData
    return (
        AnnData,
        AnnDataWrapper,
        MultivecZarrWrapper,
        VitessceConfig,
        adata_to_multivec_zarr,
        join,
        mmread,
        pd,
        vt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Load the data

        In this step, we load the raw data that has been downloaded from the HuBMAP portal https://portal.hubmapconsortium.org/browse/dataset/210d118a14c8624b6bb9610a9062656e
        """
    )
    return


@app.cell
def _(join, mmread, pd):
    mtx = mmread(join('data', 'snapatac', 'filtered_cell_by_bin.mtx')).toarray()
    barcodes_df = pd.read_csv(join('data', 'snapatac', 'barcodes.txt'), header=None)
    bins_df = pd.read_csv(join('data', 'snapatac', 'bins.txt'), header=None, names=["interval"])
    clusters_df = pd.read_csv(join('data', 'snapatac', 'umap_coords_clusters.csv'), index_col=0)
    return bins_df, clusters_df, mtx


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Convert the data to Vitessce-compatible formats

        Vitessce can load AnnData objects saved to Zarr formats efficiently.
        """
    )
    return


@app.cell
def _(bins_df):
    # The genome assembly is GRCh38 but the chromosome names in the bin names do not start with the "chr" prefix.
    # This is incompatible with the chromosome names from `negspy`, so we need to append the prefix.
    bins_df["interval"] = bins_df["interval"].apply(lambda x: "chr" + x)
    return


@app.cell
def _(AnnData, bins_df, clusters_df, mtx):
    obs = clusters_df[["cluster"]]
    obs["cluster"] = obs["cluster"].astype(str)
    obsm = { "X_umap": clusters_df[["umap.1", "umap.2"]].values }
    adata = AnnData(X=mtx, obs=obs, var=bins_df, obsm=obsm)
    adata
    return adata, obs


@app.cell
def _(join):
    multivec_zarr_path = join("data", "HBM485.TBWH.322.multivec.zarr")
    adata_zarr_path = join("data", "HBM485.TBWH.322.adata.zarr")
    return adata_zarr_path, multivec_zarr_path


@app.cell
def _(adata, adata_to_multivec_zarr, multivec_zarr_path, obs):
    # Sort cluster IDs
    cluster_ids = obs["cluster"].unique().tolist()
    cluster_ids.sort(key=int)
    # Save genomic profiles to multivec-zarr format.
    adata_to_multivec_zarr(adata, multivec_zarr_path, obs_set_col="cluster", obs_set_name="Cluster", obs_set_vals=cluster_ids)
    return


@app.cell
def _(adata, adata_zarr_path):
    # Save anndata object to AnnData-Zarr format.
    adata.write_zarr(adata_zarr_path)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 4. Make a Vitessce configuration

        We need to tell Vitessce about the data that we want to load and the visualization components that we want to include in the widget.
        For this dataset, we want to add the `GENOMIC_PROFILES` component, which renders genome browser tracks with [HiGlass](http://higlass.io).
        """
    )
    return


@app.cell
def _(
    AnnDataWrapper,
    MultivecZarrWrapper,
    VitessceConfig,
    adata_zarr_path,
    multivec_zarr_path,
    vt,
):
    vc = VitessceConfig(schema_version="1.0.15", name='HuBMAP snATAC-seq')
    dataset = vc.add_dataset(name='HBM485.TBWH.322').add_object(MultivecZarrWrapper(
        zarr_path=multivec_zarr_path
    )).add_object(AnnDataWrapper(
        adata_path=adata_zarr_path,
        obs_embedding_paths=["obsm/X_umap"],
        obs_embedding_names=["UMAP"],
        obs_set_paths=["obs/cluster"],
        obs_set_names=["Cluster"],
    ))

    genomic_profiles = vc.add_view(vt.GENOMIC_PROFILES, dataset=dataset)
    scatter = vc.add_view(vt.SCATTERPLOT, dataset=dataset, mapping = "UMAP")
    cell_sets = vc.add_view(vt.OBS_SETS, dataset=dataset)

    vc.layout(genomic_profiles / (scatter | cell_sets));
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
    vw = vc.widget(height=800)
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
