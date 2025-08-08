import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Vitessce Widget Tutorial
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of Multi-Modal Imaging Data
        We visualize IMS, PAS, and AF imaging data overlaid from the Spraggins Lab of the Biomolecular Multimodal Imaging Center (BIOMC) at Vanderbilt University, uploaded to the HuBMAP data portal.
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
    )
    from os.path import join
    return CL, VitessceConfig


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Configure Vitessce
        Set up the images from the three different assays, with the `use_physical_size_scaling` set to `True` so that the IMS image scales to the other images based on their physical sizes.
        """
    )
    return


@app.cell
def _(CL, VitessceConfig):
    vc = VitessceConfig(schema_version="1.0.16", name='Spraggins Multi-Modal', description='PAS + IMS + AF From https://portal.hubmapconsortium.org/browse/collection/6a6efd0c1a2681dc7d2faab8e4ab0bca')
    dataset = vc.add_dataset(name='Spraggins').add_file(
        url='https://assets.hubmapconsortium.org/f4188a148e4c759092d19369d310883b/ometiff-pyramids/processedMicroscopy/VAN0006-LK-2-85-PAS_images/VAN0006-LK-2-85-PAS_registered.ome.tif?token=',
        file_type="image.ome-tiff",
        coordination_values={
            "fileUid": "PAS",
        },
    )

    imageScopes = vc.add_coordination_by_dict({
        "imageLayer": CL([
          {
            "fileUid": 'PAS',
            "spatialLayerOpacity": 1,
            "spatialLayerVisible": True,
            "photometricInterpretation": 'RGB',
            "imageChannel": CL([
              {
                "spatialTargetC": 0,
                "spatialChannelColor": [255, 0, 0],
                "spatialChannelVisible": True,
                "spatialChannelOpacity": 1.0,
                "spatialChannelWindow": [0, 255],
              },
              {
                "spatialTargetC": 1,
                "spatialChannelColor": [0, 255, 0],
                "spatialChannelVisible": True,
                "spatialChannelOpacity": 1.0,
                "spatialChannelWindow": [0, 255],
              },
              {
                "spatialTargetC": 2,
                "spatialChannelColor": [0, 0, 255],
                "spatialChannelVisible": True,
                "spatialChannelOpacity": 1.0,
                "spatialChannelWindow": [0, 255],
              },
            ]),
          }
        ])
    })

    metaCoordinationScope = vc.add_meta_coordination()
    metaCoordinationScope.use_coordination_by_dict(imageScopes)

    spatial = vc.add_view("spatialBeta", dataset=dataset)
    lc = vc.add_view("layerControllerBeta", dataset=dataset)

    spatial.use_meta_coordination(metaCoordinationScope)
    lc.use_meta_coordination(metaCoordinationScope)

    vc.layout(spatial | lc);
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Create the Vitessce widget
        """
    )
    return


@app.cell
def _(vc):
    vw = vc.widget(custom_js_url='http://localhost:8000/packages/main/prod/dist/index.min.js')
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
