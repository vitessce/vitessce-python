import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import anywidget
    from vitessce import (
        VitessceConfig,
        ViewType as vt,
        CoordinationType as ct,
        ImageOmeTiffWrapper,
        ImageOmeZarrWrapper,
        AnnDataWrapper,
        SpatialDataWrapper,
        CoordinationLevel as CL,
        get_initial_coordination_scope_prefix,
        hconcat,
        vconcat,
    )
    return VitessceConfig, ImageOmeTiffWrapper, get_initial_coordination_scope_prefix


@app.cell
def _(VitessceConfig, ImageOmeTiffWrapper, get_initial_coordination_scope_prefix):
    vc = VitessceConfig(
        schema_version="1.0.17",
        name='My configuration',
        description='Data from https://portal.hubmapconsortium.org/browse/collection/6a6efd0c1a2681dc7d2faab8e4ab0bca'
    )
    dataset = vc.add_dataset(name='My dataset').add_object(
        ImageOmeTiffWrapper(
            img_url="https://assets.hubmapconsortium.org/f4188a148e4c759092d19369d310883b/ometiff-pyramids/processedMicroscopy/VAN0006-LK-2-85-PAS_images/VAN0006-LK-2-85-PAS_registered.ome.tif"
        )
    )

    spatial = vc.add_view("spatialBeta", dataset=dataset)
    lc = vc.add_view("layerControllerBeta", dataset=dataset)

    vc.link_views_by_dict([spatial, lc], {
        "photometricInterpretation": "RGB"
    }, meta=True, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))

    vc.layout(spatial | lc);
    return (vc,)


@app.cell
def _(mo, vc):
    vw = mo.ui.anywidget(vc.widget())
    vw
    return (vw,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""As you interact with the widget above, watch as the Vitessce configuration below dynamically updates. For instance, as you zoom in the spatial view, check the value of `coordinationSpace.spatialZoom`""")
    return


@app.cell
def _(vw):
    vw._config
    return


@app.cell
def _(vw):
    vw._config["coordinationSpace"]["spatialZoom"]
    return

@app.cell
def _(vc):
    import json
    print(json.dumps(vc.to_dict(), indent=2))
    return

@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
