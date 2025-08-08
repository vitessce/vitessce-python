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
        # Visualization of a SpatialData object
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
    import os
    from os.path import join, isfile, isdir
    from urllib.request import urlretrieve
    import zipfile

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
        isdir,
        isfile,
        join,
        os,
        urlretrieve,
        vt,
        zipfile,
    )


@app.cell
def _(join):
    data_dir = "data"
    zip_filepath = join(data_dir, "visium.spatialdata.zarr.zip")
    spatialdata_filepath = join(data_dir, "visium.spatialdata.zarr")
    return data_dir, spatialdata_filepath, zip_filepath


@app.cell
def _(
    data_dir,
    isdir,
    isfile,
    join,
    os,
    spatialdata_filepath,
    urlretrieve,
    zip_filepath,
    zipfile,
):
    if not isdir(spatialdata_filepath):
        if not isfile(zip_filepath):
            os.makedirs(data_dir, exist_ok=True)
            urlretrieve('https://s3.embl.de/spatialdata/spatialdata-sandbox/visium_associated_xenium_io.zip', zip_filepath)
        with zipfile.ZipFile(zip_filepath,"r") as zip_ref:
            zip_ref.extractall(data_dir)
            os.rename(join(data_dir, "data.zarr"), spatialdata_filepath)
    return


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
    vt,
):
    vc = VitessceConfig(
        schema_version="1.0.18",
        name='Visium SpatialData Demo (visium_associated_xenium_io)',
    )
    # Add data to the configuration:
    wrapper = SpatialDataWrapper(
        sdata_path=spatialdata_filepath,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        image_path="images/CytAssist_FFPE_Human_Breast_Cancer_full_image",
        table_path="tables/table",
        obs_feature_matrix_path="tables/table/X",
        obs_spots_path="shapes/CytAssist_FFPE_Human_Breast_Cancer",
        region="CytAssist_FFPE_Human_Breast_Cancer",
        coordinate_system="global",
        coordination_values={
            # The following tells Vitessce to consider each observation as a "spot"
            "obsType": "spot",
        }
    )
    dataset = vc.add_dataset(name='Breast Cancer Visium').add_object(wrapper)

    # Add views (visualizations) to the configuration:
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    feature_list = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)
    vc.link_views_by_dict([spatial, layer_controller], {
        'imageLayer': CL([{
            'photometricInterpretation': 'RGB',
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))
    obs_sets = vc.add_view(vt.OBS_SETS, dataset=dataset)
    vc.link_views([spatial, layer_controller, feature_list, obs_sets], ['obsType'], [wrapper.obs_type_label])

    # Layout the views
    vc.layout(spatial | (feature_list / layer_controller / obs_sets));
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
