import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of a SpatialData object, blobs example
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Import dependencies

        """
    )
    return


@app.cell
def _():
    import spatialdata
    from os.path import join
    return join, spatialdata


@app.cell
def _(spatialdata):
    sdata = spatialdata.datasets.blobs()
    return (sdata,)


@app.cell
def _(join):
    spatialdata_filepath = join("data", "blobs.spatialdata.zarr")
    return (spatialdata_filepath,)


@app.cell
def _(sdata, spatialdata_filepath):
    sdata.write(spatialdata_filepath, overwrite=True)
    return


@app.cell
def _():
    from vitessce import (
        VitessceConfig,
        ViewType as vt,
        CoordinationType as ct,
        CoordinationLevel as CL,
        SpatialDataWrapper,
        get_initial_coordination_scope_prefix
    )
    return (
        CL,
        SpatialDataWrapper,
        VitessceConfig,
        get_initial_coordination_scope_prefix,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Configure Vitessce

        Vitessce needs to know which pieces of data we are interested in visualizing, the visualization types we would like to use, and how we want to coordinate (or link) the views.
        """
    )
    return


@app.cell
def _(
    CL,
    SpatialDataWrapper,
    VitessceConfig,
    get_initial_coordination_scope_prefix,
    spatialdata_filepath,
):
    vc = VitessceConfig(
        schema_version="1.0.18",
        name='Visium SpatialData Demo (blobs)',
    )
    # Add data to the configuration:
    wrapper = SpatialDataWrapper(
        sdata_store=spatialdata_filepath,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        image_path="images/blobs_image",
        obs_segmentations_path="labels/blobs_labels",
        coordinate_system="global",
        coordination_values={
            "obsType": "blob",
            "fileUid": "my_unique_id"
        }
    )
    dataset = vc.add_dataset(name='Blobs').add_object(wrapper)

    # Add views (visualizations) to the configuration:
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)

    vc.link_views_by_dict([spatial, layer_controller], {
        'imageLayer': CL([{
            "fileUid": "my_unique_id",
            'photometricInterpretation': 'BlackIsZero',
            'spatialLayerOpacity': 0.9,
            'imageChannel': CL([
                {
                    "spatialTargetC": 0,
                    "spatialChannelColor": [255, 0, 0],
                    "spatialChannelOpacity": 1.0
                },
                {
                    "spatialTargetC": 1,
                    "spatialChannelColor": [0, 255, 0],
                    "spatialChannelOpacity": 1.0
                },
                {
                    "spatialTargetC": 2,
                    "spatialChannelColor": [0, 0, 255],
                    "spatialChannelOpacity": 1.0
                }
            ])
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))

    vc.link_views_by_dict([spatial, layer_controller], {
        'segmentationLayer': CL([{
            "fileUid": "my_unique_id",
            'segmentationChannel': CL([{
                'spatialChannelVisible': True,
                'obsType': 'blob',
            }]),
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"))

    # Layout the views
    vc.layout(spatial | layer_controller);
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Render the widget
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
