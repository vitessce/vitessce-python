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
        # Visualization of ComparativeData object
        """
    )
    return


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
    from oxc_py import transform
    from vitessce import VitessceConfig, hconcat, vconcat
    return VitessceConfig, hconcat, transform, vconcat


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Configure the data and views
        """
    )
    return


@app.cell
def _():
    # Reference: https://github.com/vitessce/vitessce/blob/main/examples/configs/src/view-configs/kpmp-premiere.js
    return


@app.cell
def _():
    base_url = 'https://storage.googleapis.com/vitessce-demo-data/kpmp-jan-2025/kpmp_premiere_20250330.adata.zarr'
    return (base_url,)


@app.cell
def _(VitessceConfig, base_url, hconcat, vconcat):
    vc = VitessceConfig(schema_version="1.0.17", name='Lake et al.')

    dataset = vc.add_dataset('lake_et_al').add_file(
        file_type='comparisonMetadata.anndata.zarr',
        url=base_url,
        options={
          "path": 'uns/comparison_metadata',
        },
        coordination_values={
          "obsType": 'cell',
          "sampleType": 'sample',
        },
    ).add_file(
        file_type='comparativeFeatureStats.anndata.zarr',
        url=base_url,
        options= {
          "metadataPath": 'uns/comparison_metadata',
          "indexColumn": 'names',
          "pValueColumn": 'pvals_adj',
          "foldChangeColumn": 'logfoldchanges',
          "pValueAdjusted": True,
          "foldChangeTransformation": 'log2',
        },
        coordination_values={
          "obsType": 'cell',
          "sampleType": 'sample',
          "featureType": 'gene',
        },
    ).add_file(
        file_type= 'comparativeObsSetStats.anndata.zarr',
        url= base_url,
        options= {
          "metadataPath": 'uns/comparison_metadata',
          "indexColumn": 'Cell Type',
          "interceptExpectedSampleColumn": 'Expected Sample_intercept',
          "effectExpectedSampleColumn": 'Expected Sample_effect',
          "foldChangeColumn": 'log2-fold change',
          "foldChangeTransformation": 'log2',
          "isCredibleEffectColumn": 'is_credible_effect',
        },
        coordination_values= {
          "obsType": 'cell',
          "sampleType": 'sample',
        },
    ).add_file(
      file_type='comparativeFeatureSetStats.anndata.zarr',
      url=base_url,
      options= {
        "metadataPath": 'uns/comparison_metadata',
        "indexColumn": 'pathway_name',
        "termColumn": 'pathway_term',
        "pValueColumn": 'pvals_adj',
        "pValueAdjusted": True,
        "analysisType": 'pertpy_hypergeometric',
        "featureSetLibrary": 'Reactome_2022',
      },
      coordination_values= {
        "obsType": 'cell',
        "featureType": 'gene',
        "sampleType": 'sample',
      },
    ).add_file(
      file_type='anndata.zarr',
      url=base_url,
      coordination_values={
        "obsType": 'cell',
        "featureType": 'gene',
        "featureValueType": 'expression',
        "sampleType": 'sample',
      },
      options={
        "obsFeatureMatrix": {
          "path": 'layers/pearson_residuals',
        },
        "obsEmbedding": [
          {
            "path": 'obsm/X_densmap',
            "embeddingType": 'densMAP',
          },
        ],
        "obsSets": [
          {
            "name": 'Cell Type',
            "path": 'obs/cell_type',
          },
          {
            "name": 'Subclass L1',
            "path": 'obs/subclass_l1',
          },
          {
            "name": 'Subclass L2',
            "path": 'obs/subclass_l2',
          },
          {
            "name": 'Donor ID',
            "path": 'obs/donor_id',
          },
        ],
        "sampleEdges": {
          "path": 'obs/SampleID',
        },
      },
    ).add_file(
      file_type='sampleSets.anndata.zarr',
      url=f"{base_url}/uns/__all__.samples",
      options={
        "sampleSets": [
          {
            "name": 'Disease Type',
            "path": 'diseasetype',
          },
          {
            "name": 'Adjudicated Category',
            "path": 'AdjudicatedCategory',
          },
          {
            "name": 'Enrollment Category',
            "path": 'EnrollmentCategory',
          },
        ],
      },
      coordination_values= {
        "sampleType": 'sample',
      },
    )

    biomarkerSelect = vc.add_view('biomarkerSelect', dataset=dataset, uid='biomarker-select')
    comparativeHeading = vc.add_view('comparativeHeading', dataset=dataset, uid='comparative-heading')
    dualScatterplot = vc.add_view('dualScatterplot', dataset=dataset, uid='scatterplot')
    obsSets = vc.add_view('obsSets', dataset=dataset, uid='cell-sets')
    sampleSets = vc.add_view('sampleSetPairManager', dataset=dataset, uid='sample-sets')
    obsSetSizes = vc.add_view('obsSetSizes', dataset=dataset)
    featureList = vc.add_view('featureList', dataset=dataset)
    violinPlots = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid='violin-plot')
    dotPlot = vc.add_view('dotPlot', dataset=dataset, uid='dot-plot')
    treemap = vc.add_view('treemap', dataset=dataset, uid='treemap')
    volcanoPlot = vc.add_view('volcanoPlot', dataset=dataset, uid='volcano-plot')
    volcanoPlotTable = vc.add_view('featureStatsTable', dataset=dataset, uid='volcano-plot-table')
    obsSetCompositionBarPlot = vc.add_view('obsSetCompositionBarPlot', dataset=dataset, uid='sccoda-plot')
    featureSetEnrichmentBarPlot = vc.add_view('featureSetEnrichmentBarPlot', dataset=dataset, uid='pathways-plot')

    [sampleSetScope_caseControl] = vc.add_coordination('sampleSetSelection')
    sampleSetScope_caseControl.set_value([['Disease Type', 'CKD'], ['Disease Type', 'Reference']])

    [featureSelectionScope] = vc.add_coordination('featureSelection')
    featureSelectionScope.set_value(['UMOD', 'NPHS2'])

    vc.link_views_by_dict([dualScatterplot], {
        "embeddingType": 'densMAP',
        "embeddingContoursVisible": True,
        "embeddingPointsVisible": False,
        "embeddingObsSetLabelsVisible": True,
    }, meta=False);


    vc.link_views([biomarkerSelect, dualScatterplot, obsSets, obsSetSizes, featureList, violinPlots, dotPlot, treemap, volcanoPlot, volcanoPlotTable, comparativeHeading, obsSetCompositionBarPlot, featureSetEnrichmentBarPlot, sampleSets], ['sampleType'], ['sample'])

    vc.link_views_by_dict([biomarkerSelect, dualScatterplot, obsSets, obsSetSizes, featureList, violinPlots, dotPlot, treemap, volcanoPlot, volcanoPlotTable, comparativeHeading, obsSetCompositionBarPlot, featureSetEnrichmentBarPlot, sampleSets], {
        "sampleSetSelection": sampleSetScope_caseControl,
        "featureSelection": featureSelectionScope,
    }, meta=False)

    vc.link_views_by_dict([dualScatterplot, violinPlots, featureList, dotPlot], {
        # "featureSelection": ['UMOD', 'NPHS2'], // , 'ENSG00000074803', 'ENSG00000164825'],
        "obsColorEncoding": 'geneSelection',
        "featureValueColormap": 'jet',
        "featureValueColormapRange": [0, 0.25],
        "featureAggregationStrategy": None,
    }, meta=False)

    vc.layout(hconcat(
        vconcat(dualScatterplot, biomarkerSelect, comparativeHeading, obsSets, obsSetSizes, featureList),
        vconcat(treemap, featureSetEnrichmentBarPlot, violinPlots, dotPlot, obsSetCompositionBarPlot, sampleSets),
        volcanoPlotTable,
    ));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Define the page layout
        """
    )
    return


@app.cell
def _(transform):
    PAGE_ESM = transform("""
    import clsx from "https://unpkg.com/clsx@1.1.1/dist/clsx.m.js";

    function createPage(utilsForPages) {
      const {
        React,
        usePageModeView,
      } = utilsForPages;
      function PageComponent(props) {
        const BiomarkerSelect = usePageModeView('biomarker-select');
        const ComparativeHeading = usePageModeView('comparative-heading');
        const CellSets = usePageModeView('cell-sets');
        const SampleSets = usePageModeView('sample-sets');
        const DualScatterplot = usePageModeView('scatterplot');
        const ViolinPlot = usePageModeView('violin-plot');
        const DotPlot = usePageModeView('dot-plot');
        const Treemap = usePageModeView('treemap');
        const VolcanoPlot = usePageModeView('volcano-plot');
        const VolcanoPlotTable = usePageModeView('volcano-plot-table');
        const SccodaPlot = usePageModeView('sccoda-plot');
        const PathwaysPlot = usePageModeView('pathways-plot');

    
        return (
            <>
              <style>{`
              h1, h2, h3, h4, h5, h6 {
                font-family: sans-serif;
              }
              h1 {
                font-weight: normal;
              }
              h2 {
                font-size: 36px;
              }
              h3 {
                font-size: 28px;
              }
              .stuck-comparative-heading {
                background-color: rgba(255, 255, 255, 0.7);
              }
              .stuck-comparative-heading h2 {
                font-size: 16px;
              }
              .stuck-comparative-heading h3 {
                font-size: 14px;
              }
              .view-row {
                width: 100%;
                display: flex;
                flex-direction: row;
              }
              .view-row-short {
                height: 300px;
              }
              .view-row-tall {
                height: 500px;
              }
              .view-row-left {
                width: ${(15 / 85) * 100}%;
                padding: 10px;
              }
              .view-row-left p {
                font-size: 12px;
                margin-top: 20px;
              }
              .view-row-center {
                width: ${(70 / 85) * 100}%;
              }
              .view-row-right > div {
                max-height: 50vh;
              }
              `}
              </style>
              <div style={{ width: '100%' }}>
                <div style={{ width: '70%', marginLeft: '15%' }}>
                  <h1>Comparative visualization of single-cell atlas data</h1>
                  <BiomarkerSelect />
                </div>
              </div>

              <div style={{ width: '100%', display: 'flex', flexDirection: 'row' }}>
                <div style={{ width: '85%' }}>
                  <div style={{ width: `${(70 / 85) * 100}%`, marginLeft: `${(15 / 85) * 100}%` }}>
                      <ComparativeHeading />
                  </div>
                  <div className={clsx('view-row', 'view-row-short')}>
                    <div className="view-row-left">
                      <p>This view contains a treemap visualization to communicate cell type composition in each of the selected sample groups.</p>
                    </div>
                    <div className="view-row-center">
                      <Treemap />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This view displays the results of a cell type composition analysis performed using the ScCODA algorithm (Büttner et al. 2021). Cell types with significantly different composition between the selected sample groups are displayed opaque while not-signficant results are displayed with transparent bars. The single outlined bar denotes the automatically-selected reference cell type.</p>
                    </div>
                    <div className="view-row-center">
                      <SccodaPlot />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This view displays differential expression test results, performed using the rank_genes_groups function from Scanpy (Wolf et al. 2018) with method &quot;wilcoxon&quot;. <br /><br />The arrows on the bottom left and bottom right denote the direction of the effect. Click a point in the plot to select the corresponding gene. <br /><br />Note that differential expression tests have been run for each cell type separately, so the each gene can appear multiple times (once per cell type). If there are too many points on the plot, cell types can be selected to filter the points.</p>
                    </div>
                    <div className="view-row-center">
                      <VolcanoPlot />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This view displays differential expression test results in tabular form. Click a row in the table to select the corresponding gene.</p>
                    </div>
                    <div className="view-row-center">
                      <VolcanoPlotTable />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')} style={{ height: '700px' }}>
                    <div className="view-row-left">
                      <p>This view displays gene set enrichment test results based on the differential expression results. Gene set enrichment tests have been performed using Reactome 2022 pathway gene sets from BlitzGSEA (Lachmann et al. 2022) via the hypergeometric function of Pertpy (Heumos et al. 2024).</p>
                    </div>
                    <div className="view-row-center">
                      <PathwaysPlot />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This view contains contour scatterplots which display the results of a density-preserving dimensionality reduction (Narayan et al. 2021). Contour opacities correspond to the shown percentile thresholds.</p>
                    </div>
                    <div className="view-row-center">
                      <DualScatterplot />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This dot plot view displays gene expression values per cell type and sample group for the selected biomarkers.</p>
                    </div>
                    <div className="view-row-center">
                      <DotPlot />
                    </div>
                  </div>
                  <div className={clsx('view-row', 'view-row-tall')}>
                    <div className="view-row-left">
                      <p>This violin plot view displays gene expression values per cell type and sample group for the selected biomarker.</p>
                    </div>
                    <div className="view-row-center">
                      <ViolinPlot />
                    </div>
                  </div>
                  {/* <h3>Neighborhood-level representations</h3>
                  <h1>TODO</h1>
                  <h3>Segmented instance-level representations</h3>
                  <h1>TODO</h1>
                  <h3>Image-level representations</h3>
                  <h1>TODO</h1>
                  <h3>Participant-level representations</h3>
                  <h1>TODO</h1> */}
                </div>
                <div style={{ width: '14%', marginTop: '114px', marginBottom: '100px' }}>
                    <div className="view-row-right">
                      <CellSets />
                    </div>
                    <div className="view-row-right">
                      <SampleSets />
                    </div>
                </div>

              </div>
            </>
          );
      }
      return PageComponent;
    }
    export default { createPage };
    """)
    return (PAGE_ESM,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Render page as widget
        """
    )
    return


@app.cell
def _(PAGE_ESM, vc):
    vw = vc.widget(page_esm=PAGE_ESM, page_mode=True, height=4700, prevent_scroll=False)
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
