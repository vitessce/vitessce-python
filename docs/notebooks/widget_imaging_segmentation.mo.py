import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Visualization of Segmentation Bitmask
        We visualize raw imaging data + a segmentation bitmask the [MCMicro piplene](https://mcmicro.org/) - see https://www.biorxiv.org/content/10.1101/2021.03.15.435473v1.full and specifically [Figure S1](https://www.google.com/url?q=https://www.biorxiv.org/content/10.1101/2021.03.15.435473v1.full%23F3&sa=D&source=editors&ust=1623173627976000&usg=AOvVaw3JkzCxYyE86q8jxfNCgShh)
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
    )
    from os.path import join
    return MultiImageWrapper, OmeTiffWrapper, VitessceConfig, cm


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Configure Vitessce
        Set up the two images, already pyramidal from the [bioformats2raw + raw2ometiff pipeline](https://github.com/hms-dbmi/viv/tree/master/tutorial), labeling the segmentation "on top" as the bitmask and the other as simply the image data. 
        """
    )
    return


@app.cell
def _(MultiImageWrapper, OmeTiffWrapper, VitessceConfig, cm):
    vc = VitessceConfig(schema_version="1.0.15", name='MCMicro Bitmask Visualization', description='Segmentation + Data of Exemplar 001')
    dataset = vc.add_dataset(name='MCMicro').add_object(
        MultiImageWrapper(
            image_wrappers=[
                OmeTiffWrapper(img_url='https://vitessce-demo-data.storage.googleapis.com/exemplar-001/exemplar-001.pyramid.ome.tif', name='Image'),
                OmeTiffWrapper(img_url='https://vitessce-demo-data.storage.googleapis.com/exemplar-001/cellMask.pyramid.ome.tif', name='Mask', is_bitmask=True),
            ]
     )
    )
    spatial = vc.add_view(cm.SPATIAL, dataset=dataset)
    status = vc.add_view(cm.STATUS, dataset=dataset)
    lc = vc.add_view(cm.LAYER_CONTROLLER, dataset=dataset)
    vc.layout(spatial | (lc / status));
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
    vw = vc.widget()
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
