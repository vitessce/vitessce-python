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
    sdata_url = "https://data-2.vitessce.io/data/moffitt/merfish_mouse_ileum.sdata.zarr"
    return (sdata_url,)


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
    sdata_url,
):
    vc = VitessceConfig(
        schema_version="1.0.18",
        name='SpatialData with MERFISH data',
    )
    # Add data to the configuration:

    dataset = vc.add_dataset(name='MERFISH').add_object(SpatialDataWrapper(
        sdata_url=sdata_url,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        image_path="images/stains",
        coordinate_system="global",
        coordination_values={
            "fileUid": "stains"
        }
    )).add_object(SpatialDataWrapper(
        sdata_url=sdata_url,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        obs_segmentations_path="labels/dapi_labels",
        coordinate_system="global",
        coordination_values={
            "obsType": "nucleus",
            "fileUid": "dapi"
        }
    )).add_object(SpatialDataWrapper(
        sdata_url=sdata_url,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        obs_segmentations_path="labels/membrane_labels",
        coordinate_system="global",
        coordination_values={
            "obsType": "cell",
            "fileUid": "membrane"
        }
    ))

    # Add views (visualizations) to the configuration:
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)

    vc.link_views_by_dict([spatial, layer_controller], {
        'imageLayer': CL([{
            "fileUid": "stains",
            'photometricInterpretation': 'BlackIsZero',
            'spatialLayerOpacity': 1.0,
            'spatialLayerVisible': True,
            'imageChannel': CL([
                {
                    'spatialChannelVisible': True,
                    "spatialTargetC": 0, # DAPI, Nucleus
                    "spatialChannelColor": [0, 0, 255],
                    "spatialChannelOpacity": 1.0
                },
                {
                    'spatialChannelVisible': True,
                    "spatialTargetC": 1, # Membrane, Cell
                    "spatialChannelColor": [255, 255, 255],
                    "spatialChannelOpacity": 1.0
                }
            ])
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))

    vc.link_views_by_dict([spatial, layer_controller], {
        'segmentationLayer': CL([{
            "fileUid": "membrane",
            'spatialLayerOpacity': 1.0,
            'spatialLayerVisible': True,
            'segmentationChannel': CL([{
                'spatialChannelVisible': True,
                'obsType': 'cell',
                "spatialChannelColor": [200, 200, 200],
                "obsColorEncoding": "spatialChannelColor",
            }]),
        }, {
            "fileUid": "dapi",
            'spatialLayerOpacity': 1.0,
            'spatialLayerVisible': True,
            'segmentationChannel': CL([{
                'spatialChannelVisible': True,
                'obsType': 'nucleus',
                "spatialChannelColor": [255, 255, 255],
                "obsColorEncoding": "spatialChannelColor",
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
