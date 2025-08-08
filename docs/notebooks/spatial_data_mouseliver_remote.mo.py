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
        # Visualization of a SpatialData object and individual Spatial Elements, remote data
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        This notebook explains how to create interactive visualizations of data that are accessible remotely.


        We progress through different visualization tasks, first demonstrating how Vitessce facilitates integrated imaging and spatial single-cell visualizations, then demonstrating visualization of non-spatial and image-only datasets.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        <!--## Table of Contents

        - SpatialData (via Zarr)
        - AnnData via Zarr
        - AnnData via H5AD
        - OME-Zarr
        - OME-TIFF
        -->
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Utility dependencies
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        First, we import utility dependencies which will be used to download the example dataset and manipulate file paths, zip files, and JSON files.
        """
    )
    return


@app.cell
def _():
    import os
    from os.path import join, isfile, isdir
    from urllib.request import urlretrieve
    import zipfile
    import json
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Dependencies for Vitessce
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Here, we import classes and functions from the `vitessce` Python package.
        This package includes not only APIs for [visualization configuration](https://python-docs.vitessce.io/api_config.html) but also [helper functions](https://python-docs.vitessce.io/api_data.html#vitessce-data-utils) for basic data transformation tasks.
        To specify mappings between data fields and visualization properties, the package contains [classes](https://python-docs.vitessce.io/api_data.html#module-vitessce.wrappers) which wrap standard single-cell data structures stored in formats including [AnnData](https://doi.org/10.1101/2021.12.16.473007), [SpatialData](https://doi.org/10.1038/s41592-024-02212-x), [OME-TIFF](https://doi.org/10.1007/978-3-030-23937-4_1), and [OME-Zarr](https://doi.org/10.1038/s41592-021-01326-w):

        - [AnnDataWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.AnnDataWrapper)
        - [ImageOmeTiffWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.ImageOmeTiffWrapper)
        - [ImageOmeZarrWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.ImageOmeZarrWrapper)
        - [ObsSegmentationsOmeTiffWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.ObsSegmentationsOmeTiffWrapper)
        - [ObsSegmentationsOmeZarrWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.ObsSegmentationsOmeZarrWrapper)
        - [SpatialDataWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.SpatialDataWrapper)
        """
    )
    return


@app.cell
def _():
    from vitessce import (
        VitessceConfig,
        ViewType as vt,
        CoordinationType as ct,
        CoordinationLevel as CL,
        SpatialDataWrapper,
        AnnDataWrapper,
        ImageOmeTiffWrapper,
        ImageOmeZarrWrapper,
        ObsSegmentationsOmeZarrWrapper,
        get_initial_coordination_scope_prefix,
        hconcat,
        vconcat,
    )
    from vitessce.data_utils import (
        VAR_CHUNK_SIZE,
        generate_h5ad_ref_spec,
        multiplex_img_to_ome_tiff,
        multiplex_img_to_ome_zarr,
    )
    return (
        AnnDataWrapper,
        CL,
        ImageOmeTiffWrapper,
        ImageOmeZarrWrapper,
        ObsSegmentationsOmeZarrWrapper,
        SpatialDataWrapper,
        VitessceConfig,
        get_initial_coordination_scope_prefix,
        hconcat,
        vt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        In this blog post, we perform basic data transformation tasks to save individual elements of an integrated SpatialData object to separate files in AnnData, OME-TIFF, and OME-Zarr formats.
        To perform these data transformations, we import the following dependencies.
        In general, you will typically not need to import all of these dependencies, either because you are only working with data in one of these formats, or because the data you intend to visualize is already saved to a file or directory.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Download example dataset
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        We download a mouse liver dataset which serves as a SpatialData [example dataset](https://github.com/scverse/spatialdata-notebooks/blob/main/notebooks/examples/transformations.ipynb).

        This dataset was generated by [Guilliams et al.](https://doi.org/10.1016/j.cell.2021.12.018) and processed using [SPArrOW](https://doi.org/10.1101/2024.07.04.601829) during the SpatialData [developer workshop](https://doi.org/10.37044/osf.io/8ck3e) in 2024.
        """
    )
    return


@app.cell
def _():
    base_url = "https://storage.googleapis.com/vitessce-demo-data/spatialdata-may-2025"
    return (base_url,)


@app.cell
def _(base_url):
    zip_filepath = f"{base_url}/mouse_liver.spatialdata.zarr.zip"
    spatialdata_filepath = f"{base_url}/mouse_liver.spatialdata.zarr"
    adata_zarr_filepath = f"{base_url}/mouse_liver.anndata.zarr"
    adata_h5ad_filepath = f"{base_url}/mouse_liver.h5ad"
    ref_spec_json_filepath = f"{base_url}/mouse_liver.h5ad.ref.json"
    ome_tiff_filepath = f"{base_url}/mouse_liver.ome.tif"
    offsets_json_filepath = f"{base_url}/mouse_liver.ome.tif.offsets.json"
    ome_zarr_filepath = f"{base_url}/mouse_liver.ome.zarr"
    labels_ome_zarr_filepath = f"{base_url}/mouse_liver.labels.ome.zarr"
    return (
        adata_h5ad_filepath,
        adata_zarr_filepath,
        labels_ome_zarr_filepath,
        offsets_json_filepath,
        ome_tiff_filepath,
        ome_zarr_filepath,
        ref_spec_json_filepath,
        spatialdata_filepath,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Visualization of a SpatialData object
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        SpatialData objects are the most complex type of data structure we will work with in this blog post.
        SpatialData objects function as contains for multiple types of Spatial Elements:

        - Tables (each table is represented as an AnnData object)
        - Points (e.g., coordinates of transcripts from FISH-based experiments)
        - Shapes (vector-based shapes such as polygons and circles)
        - Labels (label images, i.e., segmentation bitmasks; each label image is stored using OME-Zarr)
        - Images (microscopy images; each image is stored using OME-Zarr)
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Configure Vitessce

        Vitessce needs to know which pieces of data we are interested in visualizing, the visualization types we would like to use, and how we want to coordinate (or link) the views.
        To visualize data stored in a SpatialData object, we use the `SpatialDataWrapper` class and specify the paths (relative to the root of the Zarr [directory store](https://zarr.readthedocs.io/en/v2.18.5/api/storage.html#zarr.storage.DirectoryStore)) to different spatial elements of interest.
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
    vc = VitessceConfig(schema_version='1.0.18', name='SpatialData Demo')
    _wrapper = SpatialDataWrapper(sdata_url=spatialdata_filepath, table_path='tables/table', image_path='images/raw_image', obs_segmentations_path='labels/segmentation_mask', obs_feature_matrix_path='tables/table/X', obs_set_paths=['tables/table/obs/annotation'], obs_set_names=['Annotation'], region='nucleus_boundaries', coordinate_system='global', coordination_values={'obsType': 'cell'})
    _dataset = vc.add_dataset(name='Mouse Liver').add_object(_wrapper)
    _spatial = vc.add_view('spatialBeta', dataset=_dataset)
    feature_list = vc.add_view('featureList', dataset=_dataset)
    _layer_controller = vc.add_view('layerControllerBeta', dataset=_dataset)
    obs_sets = vc.add_view('obsSets', dataset=_dataset)
    _heatmap = vc.add_view('heatmap', dataset=_dataset)
    [obs_color_encoding_scope] = vc.add_coordination('obsColorEncoding')
    obs_color_encoding_scope.set_value('cellSetSelection')
    vc.link_views_by_dict([_spatial, _layer_controller], {'imageLayer': CL([{'photometricInterpretation': 'BlackIsZero', 'imageChannel': CL([{'spatialTargetC': 0, 'spatialChannelColor': [255, 255, 255], 'spatialChannelWindow': [0, 4000]}])}])}, scope_prefix=get_initial_coordination_scope_prefix('A', 'image'))
    vc.link_views_by_dict([_spatial, _layer_controller], {'segmentationLayer': CL([{'segmentationChannel': CL([{'obsColorEncoding': obs_color_encoding_scope}])}])}, scope_prefix=get_initial_coordination_scope_prefix('A', 'obsSegmentations'))
    vc.link_views([_spatial, _layer_controller, feature_list, obs_sets, _heatmap], ['obsType'], [_wrapper.obs_type_label])
    vc.link_views_by_dict([feature_list, obs_sets, _heatmap], {'obsColorEncoding': obs_color_encoding_scope}, meta=False)
    vc.layout(_spatial / _heatmap | _layer_controller / (feature_list | obs_sets))
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
    _vw = vc.widget()
    _vw
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Extract AnnData object from SpatialData object
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The above example demonstrates how to visualize a spatial 'omics dataset containing not only single-cell information (e.g., a cell-by-gene expression matrix, cell type annotations) but also an image and cell segmentations.
        To demonstrate how to use Vitessce to visualize data from a (non-spatial) single-cell experiment, we will extract this information from the SpatialData object and save it to a simpler [AnnData](https://anndata.readthedocs.io/) object (ignoring the imaging and spatially-resolved elements).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        As Zarr-formatted data can be easily visualized by Vitessce, we recommend saving the AnnData object to a Zarr store using the [write_zarr](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.write_zarr.html) method.
        Optionally, the shape of array chunks (for the AnnData `X` array) can be specified as a parameter, to optimize performance based on data access patterns.
        For example, a common pattern is to visualize data across all cells for one gene.
        To support such a pattern, the chunk shape can be specified as follows, `(total number of cells, small number of genes)`, resulting in tall-and-skinny array chunks.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Alternatively, your AnnData object may already be stored using the H5AD (HDF5-based) format.
        To demonstrate this scenario, we save the object using the [write_h5ad](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.write_h5ad.html) method.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        To read H5AD-formatted data, Vitessce requires an accompanying JSON [references specification](https://fsspec.github.io/kerchunk/spec.html) file, which can be constructed using the `generate_h5ad_ref_spec` utility function.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Visualization of an AnnData object
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Zarr-based AnnData
        """
    )
    return


@app.cell
def _(AnnDataWrapper, VitessceConfig, adata_zarr_filepath, vt):
    vc_1 = VitessceConfig(schema_version='1.0.18', name='AnnData (zarr)')
    _wrapper = AnnDataWrapper(adata_url=adata_zarr_filepath, obs_feature_matrix_path='X', obs_set_paths=['obs/annotation'], obs_set_names=['Annotation'], coordination_values={'obsType': 'cell'})
    _dataset = vc_1.add_dataset(name='Mouse Liver').add_object(_wrapper)
    _heatmap = vc_1.add_view(vt.HEATMAP, dataset=_dataset)
    feature_list_1 = vc_1.add_view(vt.FEATURE_LIST, dataset=_dataset)
    obs_sets_1 = vc_1.add_view(vt.OBS_SETS, dataset=_dataset)
    _violin_plots = vc_1.add_view('obsSetFeatureValueDistribution', dataset=_dataset)
    vc_1.link_views([_heatmap, feature_list_1, obs_sets_1], ['obsType', 'featureValueColormapRange'], ['cell', [0, 0.01]])
    vc_1.layout(_heatmap / _violin_plots | feature_list_1 / obs_sets_1)
    return (vc_1,)


@app.cell
def _(vc_1):
    vc_1.widget()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### H5AD-based AnnData
        """
    )
    return


@app.cell
def _(
    AnnDataWrapper,
    VitessceConfig,
    adata_h5ad_filepath,
    ref_spec_json_filepath,
    vt,
):
    vc_2 = VitessceConfig(schema_version='1.0.18', name='AnnData (h5ad)')
    _wrapper = AnnDataWrapper(adata_url=adata_h5ad_filepath, ref_url=ref_spec_json_filepath, obs_feature_matrix_path='X', obs_set_paths=['obs/annotation'], obs_set_names=['Annotation'], coordination_values={'obsType': 'cell'})
    _dataset = vc_2.add_dataset(name='Mouse Liver').add_object(_wrapper)
    _heatmap = vc_2.add_view(vt.HEATMAP, dataset=_dataset)
    feature_list_2 = vc_2.add_view(vt.FEATURE_LIST, dataset=_dataset)
    obs_sets_2 = vc_2.add_view(vt.OBS_SETS, dataset=_dataset)
    _violin_plots = vc_2.add_view('obsSetFeatureValueDistribution', dataset=_dataset)
    vc_2.link_views([_heatmap, feature_list_2, obs_sets_2], ['obsType', 'featureValueColormapRange'], ['cell', [0, 0.01]])
    vc_2.layout(_heatmap / _violin_plots | feature_list_2 / obs_sets_2)
    return feature_list_2, obs_sets_2, vc_2


@app.cell
def _(vc_2):
    vc_2.widget()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Visualization of an image file
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### OME-Zarr image
        """
    )
    return


@app.cell
def _(
    CL,
    ImageOmeZarrWrapper,
    ObsSegmentationsOmeZarrWrapper,
    VitessceConfig,
    feature_list_2,
    get_initial_coordination_scope_prefix,
    hconcat,
    labels_ome_zarr_filepath,
    obs_sets_2,
    ome_zarr_filepath,
):
    vc_3 = VitessceConfig(schema_version='1.0.18', name='Image (ome-zarr)')
    img_wrapper = ImageOmeZarrWrapper(img_url=ome_zarr_filepath, coordination_values={'fileUid': 'image'})
    segmentations_wrapper = ObsSegmentationsOmeZarrWrapper(img_url=labels_ome_zarr_filepath, coordination_values={'fileUid': 'segmentations'})
    _dataset = vc_3.add_dataset(name='Mouse Liver').add_object(img_wrapper).add_object(segmentations_wrapper)
    _spatial = vc_3.add_view('spatialBeta', dataset=_dataset)
    _layer_controller = vc_3.add_view('layerControllerBeta', dataset=_dataset)
    vc_3.link_views_by_dict([_spatial, _layer_controller], {'imageLayer': CL([{'fileUid': 'image', 'photometricInterpretation': 'BlackIsZero', 'imageChannel': CL([{'spatialTargetC': 0, 'spatialChannelColor': [255, 255, 255], 'spatialChannelWindow': [0, 4000]}])}])}, scope_prefix=get_initial_coordination_scope_prefix('A', 'image'))
    vc_3.link_views_by_dict([_spatial, _layer_controller], {'segmentationLayer': CL([{'fileUid': 'segmentations', 'segmentationChannel': CL([{'obsColorEncoding': 'spatialChannelColor', 'spatialChannelColor': [0, 255, 0], 'spatialChannelOpacity': 0.75, 'spatialSegmentationFilled': False, 'spatialSegmentationStrokeWidth': 0.25}])}])}, scope_prefix=get_initial_coordination_scope_prefix('A', 'obsSegmentations'))
    vc_3.link_views([_spatial, _layer_controller, feature_list_2, obs_sets_2], ['obsType'], ['cell'])
    vc_3.layout(hconcat(_spatial, _layer_controller, split=(2, 1)))
    return (vc_3,)


@app.cell
def _(vc_3):
    _vw = vc_3.widget()
    _vw
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### OME-TIFF image
        """
    )
    return


@app.cell
def _(
    CL,
    ImageOmeTiffWrapper,
    VitessceConfig,
    get_initial_coordination_scope_prefix,
    hconcat,
    offsets_json_filepath,
    ome_tiff_filepath,
):
    vc_4 = VitessceConfig(schema_version='1.0.18', name='Image and segmentations (ome-tiff)')
    _wrapper = ImageOmeTiffWrapper(img_url=ome_tiff_filepath, offsets_url=offsets_json_filepath)
    _dataset = vc_4.add_dataset(name='Mouse Liver').add_object(_wrapper)
    _spatial = vc_4.add_view('spatialBeta', dataset=_dataset)
    _layer_controller = vc_4.add_view('layerControllerBeta', dataset=_dataset)
    vc_4.link_views_by_dict([_spatial, _layer_controller], {'imageLayer': CL([{'photometricInterpretation': 'BlackIsZero', 'imageChannel': CL([{'spatialTargetC': 0, 'spatialChannelColor': [255, 255, 255], 'spatialChannelWindow': [0, 4000]}])}])}, scope_prefix=get_initial_coordination_scope_prefix('A', 'image'))
    vc_4.layout(hconcat(_spatial, _layer_controller, split=(2, 1)))
    return (vc_4,)


@app.cell
def _(vc_4):
    vc_4.widget()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Data location options

        Vitessce can visualize data [not only stored](https://python-docs.vitessce.io/data_options.html) locally (referenced using local file or directory paths) but also stored remotely (referenced using absolute URL paths).
        Depending on whether the Python kernel running the Jupyter process is running locally versus remotely (e.g., on a cluster or cloud platform such as Google Colab), certain data locations may be challenging to access from the machine running the Jupyter notebook frontend (i.e., the machine on which the web browser used to view the notebook is installed).

        To provide data via one of these alternative mechanisms, use parameters with the following suffices when instantiating the data wrapper classes.

        - `_path`: Local file or directory
        - `_url`: Remote file or directory
        - `_store`: Zarr-store-accessible (for zarr-based formats)
        - `_artifact`: Lamin artifact

        For example, `adata_path` can be exchanged for one of the following options.

        ```diff
        AnnDataWrapper(
        -    adata_path="./mouse_liver.spatialdata.zarr",
        +    adata_url="https://example.com/mouse_liver.spatialdata.zarr", # Absolute URL
        +    adata_store="./mouse_liver.spatialdata.zarr", # String interpreted as root of DirectoryStore
        +    adata_store=zarr.DirectoryStore("./mouse_liver.spatialdata.zarr"), # Instance of zarr.storage
        +    adata_artifact=adata_zarr_artifact, # Instance of ln.Artifact
            ...
        ```

        Note that the `_store` options are only available for Zarr-based formats, such as AnnDataWrapper, SpatialDataWrapper, and ImageOmeZarrWrapper.
        Further, when multiple files are required, such as both `adata_path` and `ref_path` (for the JSON reference specification file accompanying an H5AD file) or `img_path` and `offsets_path`, multiple parameter suffices may need to be changed.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
