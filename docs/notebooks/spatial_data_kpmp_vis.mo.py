import marimo

__generated_with = "0.13.15"
app = marimo.App()


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
    # Reference: https://github.com/vitessce/vitessce/blob/main/examples/configs/src/view-configs/spatial-beta/kpmp.js
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
    from os.path import join
    return (
        CL,
        SpatialDataWrapper,
        VitessceConfig,
        get_initial_coordination_scope_prefix,
        hconcat,
        vconcat,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Configure Vitessce
        """
    )
    return


@app.cell
def _():
    sdata_url = "https://storage.googleapis.com/vitessce-demo-data/kpmp-f2f-march-2023/S-1905-017737/sdata.zarr"
    return (sdata_url,)


@app.cell
def _(
    CL,
    SpatialDataWrapper,
    VitessceConfig,
    get_initial_coordination_scope_prefix,
    hconcat,
    sdata_url,
    vconcat,
):
    # Create a VitessceConfig instance.
    vc = VitessceConfig(schema_version="1.0.18", name="SpatialData")

    t_obstype = "Tubule"
    a_obstype = "Artery"
    ci_obstype = "Cortical Interstitium"
    gsg_obstype = "G. S. Glomerulus"
    ngsg_obstype = "Non-G. S. Glomerulus"
    ifta_obstype = "Interstitial Fibrosis and Tubular Atrophy"
    ptc_obstype = "Peritubular Capillaries"

    # Add a new dataset to the Vitessce configuration,
    # then add the wrapper class instance to this dataset.
    dataset = vc.add_dataset(name='KPMP').add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            image_path="images/image",
            coordinate_system="global",
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_tubules",
            obs_segmentations_path="labels/labels_tubules",
            obs_feature_matrix_path="tables/table_tubules/X",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_tubules",
              "obsType": t_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            obs_segmentations_path="labels/labels_arteries_arterioles",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_arteries_arterioles",
              "obsType": a_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_cortical_interstitia",
            obs_segmentations_path="labels/labels_cortical_interstitia",
            obs_feature_matrix_path="tables/table_cortical_interstitia/X",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_cortical_interstitia",
              "obsType": ci_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_globally_sclerotic_glomeruli",
            obs_segmentations_path="labels/labels_globally_sclerotic_glomeruli",
            obs_feature_matrix_path="tables/table_globally_sclerotic_glomeruli/X",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_globally_sclerotic_glomeruli",
              "obsType": gsg_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_non_globally_sclerotic_glomeruli",
            obs_segmentations_path="labels/labels_non_globally_sclerotic_glomeruli",
            obs_feature_matrix_path="tables/table_non_globally_sclerotic_glomeruli/X",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_non_globally_sclerotic_glomeruli",
              "obsType": ngsg_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_interstitialfibrosis_and_tubular_atrophy",
            obs_segmentations_path="labels/labels_interstitialfibrosis_and_tubular_atrophy",
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_interstitialfibrosis_and_tubular_atrophy",
              "obsType": ifta_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    ).add_object(
        SpatialDataWrapper(
            sdata_url=sdata_url,
            table_path="tables/table_peritubular_capillaries",
            obs_segmentations_path="labels/labels_peritubular_capillaries",
            obs_feature_matrix_path="tables/table_peritubular_capillaries/X",
            obs_set_paths=[
                "tables/table_peritubular_capillaries/obs/cortex_ifta_set",
                "tables/table_peritubular_capillaries/obs/cortex_set",
                "tables/table_peritubular_capillaries/obs/ifta_set",
                ["tables/table_peritubular_capillaries/obs/cortex_set", "tables/table_peritubular_capillaries/obs/ifta_set"],
            ],
            obs_set_names=[
                "Cortex and IFTA membership",
                "Cortex membership",
                "IFTA membership",
                "Cortex and IFTA hierarchy",
            ],
            coordinate_system="global",
            coordination_values={
              "fileUid": "labels_peritubular_capillaries",
              "obsType": ptc_obstype,
              "featureType": 'feature',
              "featureValueType": 'value',
            }
        )
    )

    # Add views (visualizations) to the configuration.
    spatial = vc.add_view("spatialBeta", dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)
    # Tubules
    tubules_feature_list = vc.add_view("featureList", dataset=dataset).set_props(title="Tubules")
    tubules_histogram = vc.add_view("featureValueHistogram", dataset=dataset)
    # Peritubular capillaries
    pt_feature_list = vc.add_view("featureList", dataset=dataset).set_props(title="Peritubular Capillaries")
    pt_histogram = vc.add_view("featureValueHistogram", dataset=dataset)
    # GSG
    gsg_feature_list = vc.add_view("featureList", dataset=dataset).set_props(title="Globally Sclerotic Glomeruli")
    gsg_histogram = vc.add_view("featureValueHistogram", dataset=dataset)
    # NGSG
    ngsg_feature_list = vc.add_view("featureList", dataset=dataset).set_props(title="Non-Globally Sclerotic Glomeruli")
    ngsg_histogram = vc.add_view("featureValueHistogram", dataset=dataset)

    # Add obsSets, obsSetSizes, and violin plot views for PTC sets+areas/aspectRatio
    pt_sets = vc.add_view("obsSets", dataset=dataset)
    pt_bar_plot = vc.add_view("obsSetSizes", dataset=dataset)
    pt_violin_plot = vc.add_view("obsSetFeatureValueDistribution", dataset=dataset).set_props(jitter=True)


    # Coordination of views.
    [ft_scope, fvt_scope] = vc.add_coordination("featureType", "featureValueType")
    ft_scope.set_value("feature")
    fvt_scope.set_value("value")

    [t_ot_scope, t_fs_scope, t_oce_scope] = vc.add_coordination("obsType", "featureSelection", "obsColorEncoding")
    t_ot_scope.set_value(t_obstype)
    t_oce_scope.set_value("spatialChannelColor")

    [pt_ot_scope, pt_fs_scope, pt_oce_scope] = vc.add_coordination("obsType", "featureSelection", "obsColorEncoding")
    pt_ot_scope.set_value(ptc_obstype)
    pt_fs_scope.set_value(["Area"])
    pt_oce_scope.set_value("cellSetSelection")


    [gsg_ot_scope, gsg_fs_scope, gsg_oce_scope, gsg_fvcr_scope] = vc.add_coordination("obsType", "featureSelection", "obsColorEncoding", "featureValueColormapRange")
    gsg_ot_scope.set_value(gsg_obstype)
    gsg_oce_scope.set_value("spatialChannelColor")
    gsg_fvcr_scope.set_value([(3077 - 2333) / (29911 - 2333), 1.0])

    [ngsg_ot_scope, ngsg_fs_scope, ngsg_oce_scope, ngsg_fvcr_scope] = vc.add_coordination("obsType", "featureSelection", "obsColorEncoding", "featureValueColormapRange")
    ngsg_ot_scope.set_value(ngsg_obstype)
    ngsg_oce_scope.set_value("spatialChannelColor")
    ngsg_fvcr_scope.set_value([0.0, 1 - (59451 - 29911) / (59451 - 3077)])

    vc.link_views_by_dict([spatial, layer_controller], {
        "imageLayer": CL([{
            "spatialLayerOpacity": 0.1,
            "photometricInterpretation": "RGB",
        }]),
    }, meta=True, scope_prefix=get_initial_coordination_scope_prefix("A", "image"))

    vc.link_views_by_dict([spatial, layer_controller], {
        "segmentationLayer": CL([
            {
                "fileUid": "labels_tubules",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": t_ot_scope,
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "featureSelection": t_fs_scope,
                    "spatialChannelVisible": False,
                    "spatialChannelColor": [73, 155, 119],
                    "spatialChannelOpacity": 0.5,
                    "obsColorEncoding": t_oce_scope,
                    "featureValueColormapRange": [0, 1],
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_arteries_arterioles",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": "Artery",
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "spatialChannelVisible": False,
                    "spatialChannelColor": [237, 226, 107],
                    "spatialChannelOpacity": 0.5,
                    "obsColorEncoding": "spatialChannelColor",
                    "featureValueColormapRange": [0, 1],
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_cortical_interstitia",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": "Cortical Interstitium",
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "spatialChannelVisible": False,
                    "spatialChannelColor": [255, 255, 255],
                    "spatialChannelOpacity": 0.5,
                    "obsColorEncoding": "spatialChannelColor",
                    "featureValueColormapRange": [0, 1],
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_globally_sclerotic_glomeruli",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": gsg_ot_scope,
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "featureSelection": gsg_fs_scope,
                    "spatialChannelVisible": False,
                    "spatialChannelColor": [52, 113, 171],
                    "spatialChannelOpacity": 0.5,
                    "obsColorEncoding": gsg_oce_scope,
                    "featureValueColormapRange": gsg_fvcr_scope,
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_non_globally_sclerotic_glomeruli",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": ngsg_ot_scope,
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "featureSelection": ngsg_fs_scope,
                    "spatialChannelVisible": False,
                    "spatialChannelColor": [114, 179, 226],
                    "spatialChannelOpacity": 0.5,
                    "obsColorEncoding": ngsg_oce_scope,
                    "featureValueColormapRange": ngsg_fvcr_scope,
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_interstitialfibrosis_and_tubular_atrophy",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": "Interstitial Fibrosis and Tubular Atrophy",
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "spatialChannelVisible": True,
                    "spatialChannelColor": [218, 161, 66],
                    "spatialChannelOpacity": 1.0,
                    "obsColorEncoding": "spatialChannelColor",
                    "featureValueColormapRange": [0, 1],
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": False,
                    "obsHighlight": None,
                }]),
            },
            {
                "fileUid": "labels_peritubular_capillaries",
                "segmentationChannel": CL([{
                    "spatialTargetC": 0,
                    "obsType": pt_ot_scope,
                    "featureType": ft_scope,
                    "featureValueType": fvt_scope,
                    "featureSelection": pt_fs_scope,
                    "spatialChannelVisible": True,
                    "spatialChannelColor": [197, 101, 47],
                    "spatialChannelOpacity": 1.0,
                    "obsColorEncoding": pt_oce_scope,
                    "featureValueColormapRange": [0, 0.5],
                    "featureAggregationStrategy": "first",
                    "spatialSegmentationFilled": True,
                    "obsHighlight": None,
                }]),
            }
        ]),
    }, meta=True, scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"))

    tubules_feature_list.use_coordination(t_ot_scope, ft_scope, fvt_scope, t_fs_scope, t_oce_scope)
    tubules_histogram.use_coordination(t_ot_scope, ft_scope, fvt_scope, t_fs_scope, t_oce_scope)

    pt_feature_list.use_coordination(pt_ot_scope, ft_scope, fvt_scope, pt_fs_scope, pt_oce_scope)
    pt_histogram.use_coordination(pt_ot_scope, ft_scope, fvt_scope, pt_fs_scope, pt_oce_scope)

    pt_sets.use_coordination(pt_ot_scope, ft_scope, fvt_scope, pt_fs_scope, pt_oce_scope)
    pt_bar_plot.use_coordination(pt_ot_scope, ft_scope, fvt_scope, pt_fs_scope, pt_oce_scope)
    pt_violin_plot.use_coordination(pt_ot_scope, ft_scope, fvt_scope, pt_fs_scope, pt_oce_scope)


    gsg_feature_list.use_coordination(gsg_ot_scope, ft_scope, fvt_scope, gsg_fs_scope, gsg_oce_scope, gsg_fvcr_scope)
    gsg_histogram.use_coordination(gsg_ot_scope, ft_scope, fvt_scope, gsg_fs_scope, gsg_oce_scope, gsg_fvcr_scope)

    ngsg_feature_list.use_coordination(ngsg_ot_scope, ft_scope, fvt_scope, ngsg_fs_scope, ngsg_oce_scope, ngsg_fvcr_scope)
    ngsg_histogram.use_coordination(ngsg_ot_scope, ft_scope, fvt_scope, ngsg_fs_scope, ngsg_oce_scope, ngsg_fvcr_scope)


    # Layout the views in a grid arrangement.
    vc.layout(vconcat(
        hconcat(spatial, layer_controller, split=[3, 1]),
        hconcat(
            (tubules_feature_list / tubules_histogram),
            (pt_feature_list / pt_histogram),
            (gsg_feature_list / gsg_histogram),
            (ngsg_feature_list / ngsg_histogram)
        ),
        hconcat(pt_sets, pt_bar_plot, pt_violin_plot)
    ));
    return (vc,)


@app.cell
def _(vc):
    vw = vc.widget(height=1000)
    vw
    return


@app.cell
def _():
    #import json
    #print(json.dumps(vc.to_dict(), indent=2))
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
