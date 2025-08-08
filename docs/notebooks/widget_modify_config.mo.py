import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Programmatic modification of widget configuration
        """
    )
    return


@app.cell
def _():
    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        OmeTiffWrapper,
        MultiImageWrapper,
        CoordinationLevel as CL,
        ObsSegmentationsOmeTiffWrapper,
        ImageOmeTiffWrapper,
        get_initial_coordination_scope_prefix,
    )
    import random
    return (
        CL,
        ImageOmeTiffWrapper,
        ObsSegmentationsOmeTiffWrapper,
        VitessceConfig,
        get_initial_coordination_scope_prefix,
        random,
    )


@app.cell
def _(
    CL,
    ImageOmeTiffWrapper,
    ObsSegmentationsOmeTiffWrapper,
    VitessceConfig,
    get_initial_coordination_scope_prefix,
):
    vc = VitessceConfig(schema_version="1.0.16")
    dataset = vc.add_dataset(name='Spraggins').add_object(
        ImageOmeTiffWrapper(
            img_url="https://storage.googleapis.com/vitessce-demo-data/kpmp-f2f-march-2023/S-1905-017737/S-1905-017737_PAS_2of2_bf.ome.tif",
            offsets_url="https://storage.googleapis.com/vitessce-demo-data/kpmp-f2f-march-2023/S-1905-017737/S-1905-017737_PAS_2of2_bf.offsets.json"
        )
    ).add_object(
        ObsSegmentationsOmeTiffWrapper(
            img_url="https://storage.googleapis.com/vitessce-demo-data/kpmp-f2f-march-2023/S-1905-017737/S-1905-017737_PAS_2of2.ome.tif",
            offsets_url="https://storage.googleapis.com/vitessce-demo-data/kpmp-f2f-march-2023/S-1905-017737/S-1905-017737_PAS_2of2.offsets.json",
            obs_types_from_channel_names=True
        )
    )

    spatial = vc.add_view("spatialBeta", dataset=dataset)
    lc = vc.add_view("layerControllerBeta", dataset=dataset)

    vc.link_views_by_dict([spatial, lc], {
        "imageLayer": CL([
            {
              "photometricInterpretation": "RGB"
            }
        ]),
    }, meta=True, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))

    vc.layout(spatial | lc);
    return (vc,)


@app.cell
def _(vc):
    vw = vc.widget(remount_on_uid_change=False)
    vw
    return (vw,)


@app.cell
def _(vw):
    # Inspect the current configuration value.
    # This is a dict in the JSON-based format https://vitessce.io/docs/view-config-json/
    vw.config
    return


@app.cell
def _(random, vw):
    # Programatically set a different zoom level and toggle the visibility/color of different segmentation layers:
    vw.config = {
        **vw.config,
        # Need to provide a fresh "uid" value.
        # This will tell Vitessce that the contents should be diff-ed against the previous config.
        "uid": f"new_config_{random.random()}",
        "coordinationSpace": {
          # Information about the coordination space can be found at https://vitessce.io/docs/coordination-types/
          **vw.config["coordinationSpace"],
          "spatialZoom": {
              **vw.config["coordinationSpace"]["spatialZoom"],
              "A": -8
          },
          "spatialChannelVisible": {
              **vw.config["coordinationSpace"]["spatialChannelVisible"],
              "init_A_obsSegmentations_0": True,
              "init_A_obsSegmentations_1": False,
              "init_A_obsSegmentations_2": False,
              "init_A_obsSegmentations_3": False
          },
          "spatialChannelColor": {
              **vw.config["coordinationSpace"]["spatialChannelColor"],
              "init_A_obsSegmentations_0": [255, 0, 0],
          }
      }
    }
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
