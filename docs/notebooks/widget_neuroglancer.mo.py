import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Example usage of Neuroglancer view
        """
    )
    return


@app.cell
def _():
    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        ImageOmeTiffWrapper,
        CsvWrapper,
        hconcat,
        vconcat,
        get_initial_coordination_scope_prefix,
        CoordinationLevel as CL
    )
    from os.path import join
    return (
        CL,
        CsvWrapper,
        ImageOmeTiffWrapper,
        VitessceConfig,
        get_initial_coordination_scope_prefix,
        hconcat,
        vconcat,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Configure Vitessce
        """
    )
    return


@app.cell
def _(
    CL,
    CsvWrapper,
    ImageOmeTiffWrapper,
    VitessceConfig,
    get_initial_coordination_scope_prefix,
    hconcat,
    vconcat,
):
    vc = VitessceConfig(schema_version="1.0.17")
    dataset = vc.add_dataset(name='Meshes').add_object(
        ImageOmeTiffWrapper(
            img_url='https://lsp-public-data.s3.amazonaws.com/yapp-2023-3d-melanoma/Dataset1-LSP13626-invasive-margin.ome.tiff',
            offsets_url='https://lsp-public-data.s3.amazonaws.com/yapp-2023-3d-melanoma/Dataset1-LSP13626-invasive-margin.offsets.json',
            coordination_values={
              "fileUid": 'melanoma',
            },
        )
    ).add_object(
        CsvWrapper(
            data_type="obsEmbedding",
            csv_url='https://storage.googleapis.com/vitessce-demo-data/neuroglancer-march-2025/melanoma_with_embedding_filtered_ids.csv',
            options= {
              "obsIndex": 'id',
              "obsEmbedding": ['tSNE1', 'tSNE2'],
            },
            coordination_values= {
              "obsType": 'cell',
              "embeddingType": 'TSNE',
            },
        )
    ).add_object(
        CsvWrapper(
            data_type="obsSets",
            csv_url='https://storage.googleapis.com/vitessce-demo-data/neuroglancer-march-2025/melanoma_with_embedding_filtered_ids.csv',
            coordination_values={
              "obsType": 'cell',
            },
            options= {
              "obsIndex": 'id',
              "obsSets": [
                {
                  "name": 'Clusters',
                  "column": 'cluster',
                },
              ],
            },
        )
    )
    spatialThreeView = vc.add_view('spatialBeta', dataset=dataset);
    lcView = vc.add_view('layerControllerBeta', dataset=dataset);
    obsSets = vc.add_view('obsSets', dataset=dataset);
    scatterView = vc.add_view('scatterplot', dataset=dataset, mapping="TSNE");
    # Configuration via props.viewerState is temporary and subject to change.
    neuroglancerView = vc.add_view('neuroglancer', dataset=dataset).set_props(viewerState={
        "dimensions": {
          "x": [
            1e-9,
            "m"
          ],
          "y": [
            1e-9,
            "m"
          ],
          "z": [
            1e-9,
            "m"
          ]
        },
        "position": [
          49.5,
          1000.5,
          5209.5
        ],
        "crossSectionScale": 1,
        "projectionOrientation": [
          -0.636204183101654,
          -0.5028395652770996,
          0.5443811416625977,
          0.2145828753709793
        ],
        "projectionScale": 1024,
        "layers": [
          {
            "type": "segmentation",
            "source": "precomputed://https://vitessce-data-v2.s3.us-east-1.amazonaws.com/data/sorger/invasive_meshes",
            "segments": [
              "5"
            ],
            "segmentColors": {
              "5": "red"
            },
            "name": "segmentation"
          }
        ],
        "showSlices": False,
        "layout": "3d"
    });

    vc.link_views([scatterView], ['embeddingObsRadiusMode', 'embeddingObsRadius'], ['manual', 4]);

    # Sync the zoom/rotation/pan states
    vc.link_views_by_dict([spatialThreeView, lcView, neuroglancerView], {
        "spatialRenderingMode": '3D',
        "spatialZoom": 0,
        "spatialTargetT": 0,
        "spatialTargetX": 0,
        "spatialTargetY": 0,
        "spatialTargetZ": 0,
        "spatialRotationX": 0,
        "spatialRotationY": 0,
    }, meta=False);

    # Initialize the image properties
    vc.link_views_by_dict([spatialThreeView, lcView], {
        "imageLayer": CL([
          {
            "fileUid": 'melanoma',
            "spatialLayerOpacity": 1,
            "spatialTargetResolution": None,
            "imageChannel": CL([
              {
                "spatialTargetC": 0,
                "spatialChannelColor": [255, 0, 0],
                "spatialChannelVisible": True,
                "spatialChannelOpacity": 1.0,
              },
            ]),
          },
        ]),
    }, scope_prefix=get_initial_coordination_scope_prefix('A', 'image'));


    vc.layout(hconcat(neuroglancerView, spatialThreeView, vconcat(lcView, obsSets, scatterView)));
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
