import marimo

__generated_with = "0.13.15"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Vitessce Widget Tutorial""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Visualization of a SpatialData object""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Import dependencies""")
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
    zip_filepath = join(data_dir, "visium_hd_3.0.0_io.zip")
    spatialdata_filepath = join(data_dir, "visium_hd_3.0.0.spatialdata.zarr")
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
            urlretrieve('https://s3.embl.de/spatialdata/spatialdata-sandbox/visium_hd_3.0.0_io.zip', zip_filepath)
        with zipfile.ZipFile(zip_filepath,"r") as zip_ref:
            zip_ref.extractall(data_dir)
            os.rename(join(data_dir, "data.zarr"), spatialdata_filepath)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Rasterize bins
    Reference: https://spatialdata.scverse.org/en/stable/tutorials/notebooks/notebooks/examples/technology_visium_hd.html#performant-on-the-fly-data-rasterization
    """
    )
    return


@app.cell
def _():
    from spatialdata import (
        read_zarr,
        rasterize_bins,
        rasterize_bins_link_table_to_labels
    )
    return rasterize_bins, rasterize_bins_link_table_to_labels, read_zarr


@app.cell
def _(
    rasterize_bins,
    rasterize_bins_link_table_to_labels,
    read_zarr,
    spatialdata_filepath,
):
    sdata = read_zarr(spatialdata_filepath)

    for bin_size in ["016", "008", "002"]:
        # rasterize_bins() requires a compresed sparse column (csc) matrix
        sdata.tables[f"square_{bin_size}um"].X = sdata.tables[f"square_{bin_size}um"].X.tocsc()
        rasterized = rasterize_bins(
            sdata,
            f"Visium_HD_Mouse_Small_Intestine_square_{bin_size}um",
            f"square_{bin_size}um",
            "array_col",
            "array_row",
            # We want to rasterize to a Labels element, rather than an Image element.
            return_region_as_labels=True
        )
        sdata[f"rasterized_{bin_size}um"] = rasterized
        rasterize_bins_link_table_to_labels(
            sdata,
            table_name=f"square_{bin_size}um",
            rasterized_labels_name=f"rasterized_{bin_size}um",
        )
        try:
            sdata.write_element(f"rasterized_{bin_size}um")
        except:
            pass
        
    sdata
    return (sdata,)


@app.cell
def _(sdata):
    sdata
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
        name='Visium HD SpatialData Demo',
    )
    # Add data to the configuration:
    wrapper = SpatialDataWrapper(
        sdata_path=spatialdata_filepath,
        # The following paths are relative to the root of the SpatialData zarr store on-disk.
        image_path="images/Visium_HD_Mouse_Small_Intestine_full_image",
        table_path="tables/square_016um",
        obs_feature_matrix_path="tables/square_016um/X",
        obs_segmentations_path="labels/rasterized_016um",
        #region="CytAssist_FFPE_Human_Breast_Cancer",
        coordinate_system="Visium_HD_Mouse_Small_Intestine",
        coordination_values={
            # The following tells Vitessce to consider each observation as a "bin"
            "obsType": "bin",
        }
    )
    dataset = vc.add_dataset(name='Visium HD').add_object(wrapper)

    # Add views (visualizations) to the configuration:
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    feature_list = vc.add_view("featureList", dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)
    vc.link_views_by_dict([spatial, layer_controller], {
        'imageLayer': CL([{
            'photometricInterpretation': 'RGB',
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))
    vc.link_views_by_dict([spatial, layer_controller], {
        'segmentationLayer': CL([{
            'segmentationChannel': CL([{
                'spatialChannelOpacity': 0.5,
                'obsColorEncoding': 'geneSelection',
                'featureValueColormapRange': [0, 0.5],
            }])
        }]),
    }, scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"))
    obs_sets = vc.add_view(vt.OBS_SETS, dataset=dataset)
    vc.link_views([spatial, layer_controller, feature_list, obs_sets], ['obsType', 'featureSelection'], [wrapper.obs_type_label, ['AA986860']])

    # Layout the views
    vc.layout(spatial | (feature_list / layer_controller / obs_sets));
    return (vc,)


@app.cell
def _(mo):
    mo.md(r"""### Render the widget""")
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
