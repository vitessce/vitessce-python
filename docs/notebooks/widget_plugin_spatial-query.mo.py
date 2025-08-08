import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Vitessce SpatialQuery plugin usage demo
        """
    )
    return


app._unparsable_cell(
    r"""
    #!pip install \"vitessce[all]==3.3.0\" esbuild_py anndata
    !pip install \"mlxtend~=0.23.0\"
    #!pip install -i \"https://test.pypi.org/simple/\" SpatialQuery
    !pip install \"SpatialQuery @ git+https://github.com/ShaokunAn/Spatial-Query@main\"
    """,
    name="_"
)


@app.cell
def _():
    from os.path import join
    from anndata import read_h5ad
    from vitessce import (
        VitessceConfig,
        AnnDataWrapper,
        ViewType as vt,
        CoordinationType as ct,
        CoordinationLevel as CL,
    )
    from vitessce.widget_plugins import SpatialQueryPlugin
    return (
        AnnDataWrapper,
        CL,
        SpatialQueryPlugin,
        VitessceConfig,
        join,
        read_h5ad,
    )


@app.cell
def _(join, read_h5ad):
    adata = read_h5ad(join("data", "HBM987_KWLK_254", "secondary_analysis.h5ad"))
    zarr_path = join("data", "HBM987_KWLK_254", "secondary_analysis.h5ad.zarr")
    adata.write_zarr(zarr_path)
    return adata, zarr_path


@app.cell
def _(SpatialQueryPlugin, adata):
    plugin = SpatialQueryPlugin(adata)
    return (plugin,)


@app.cell
def _(AnnDataWrapper, CL, VitessceConfig, plugin, zarr_path):
    vc = VitessceConfig(schema_version="1.0.16", name="Spatial-Query")
    dataset = vc.add_dataset("Query results").add_object(AnnDataWrapper(
        adata_path=zarr_path,
        obs_feature_matrix_path="X",
        obs_set_paths=["obs/predicted.ASCT.celltype"],
        obs_set_names=["Cell Type"],
        obs_spots_path="obsm/X_spatial",
        feature_labels_path="var/hugo_symbol",
        coordination_values={
            "featureLabelsType": "Gene symbol",
        }
    ))

    spatial_view = vc.add_view("spatialBeta", dataset=dataset)
    lc_view = vc.add_view("layerControllerBeta", dataset=dataset)
    sets_view = vc.add_view("obsSets", dataset=dataset)
    features_view = vc.add_view("featureList", dataset=dataset)
    sq_view = vc.add_view("spatialQuery", dataset=dataset)

    obs_set_selection_scope, = vc.add_coordination("obsSetSelection",)
    obs_set_selection_scope.set_value(None)

    sets_view.use_coordination(obs_set_selection_scope)
    sq_view.use_coordination(obs_set_selection_scope)
    spatial_view.use_coordination(obs_set_selection_scope)
    features_view.use_coordination(obs_set_selection_scope)

    vc.link_views([spatial_view, lc_view, sets_view, features_view],
        ["additionalObsSets", "obsSetColor"],
        [plugin.additional_obs_sets, plugin.obs_set_color]
    )
    vc.link_views_by_dict([spatial_view, lc_view], {
        "spotLayer": CL([
            {
                "obsType": "cell",
                "spatialSpotRadius": 15,
            },
        ])
    })

    vc.layout((spatial_view | (lc_view / features_view)) / (sets_view | sq_view));
    return (vc,)


@app.cell
def _(plugin, vc):
    vw = vc.widget(height=900, plugins=[plugin], remount_on_uid_change=False)
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
