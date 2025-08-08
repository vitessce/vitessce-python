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
        # Page mode example
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
    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        AnnDataWrapper,
        CsvWrapper,
    )
    from oxc_py import transform
    return AnnDataWrapper, VitessceConfig, cm, transform


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
    url = 'https://storage.googleapis.com/vitessce-demo-data/anndata-test/pbmc3k_processed.zarr'
    return (url,)


@app.cell
def _(AnnDataWrapper, VitessceConfig, cm, url):
    vc = VitessceConfig(schema_version="1.0.17", name='PBMC Reference')
    dataset = vc.add_dataset(name='PBMC 3k').add_object(
        AnnDataWrapper(
            adata_url=url,
            obs_set_paths=["obs/louvain"],
            obs_set_names=["Louvain"],
            obs_embedding_paths=["obsm/X_umap", "obsm/X_pca"],
            obs_embedding_names=["UMAP", "PCA"],
            obs_feature_matrix_path="X"
        )
    )

    umap = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP", uid="scatterplot-umap")
    pca = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="PCA", uid="scatterplot-pca")
    cell_sets = vc.add_view(cm.OBS_SETS, dataset=dataset, uid="cell-sets")
    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset, uid="gene-list")
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset, uid="heatmap")

    vc.layout((umap / pca) | ((cell_sets | genes) / heatmap));
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
    function createPage(utilsForPages) {
      const {
        React,
        usePageModeView,
      } = utilsForPages;
      function PageComponent(props) {
        const ScatterplotUmap = usePageModeView('scatterplot-umap');
        const ScatterplotPca = usePageModeView('scatterplot-pca');
        const CellSets = usePageModeView('cell-sets');
        const GeneList = usePageModeView('gene-list');
        const Heatmap = usePageModeView('heatmap');
    
        return (
            <>
              <style>{`
              h1, h2, h3, h4, h5, h6 {
                font-family: sans-serif;
              }
              h3 {
                font-size: 20px;
              }
              .fancy-heading {
                  text-shadow: 1px 1px 2px pink;
              }
              `}
              </style>
              <div style={{ width: '100%', display: 'flex', flexDirection: 'row', background: 'lightblue' }}>
                <div style={{ width: '80%'}}>
                  <h3 style={{ fontFamily: 'Courier New' }} className="fancy-heading">This is an arbitrary HTML element with custom CSS</h3>
                  <div style={{ width: '100%', height: '400px', display: 'flex', flexDirection: 'row' }}>
                    <div style={{ width: '45%' }}>
                      <ScatterplotUmap />
                    </div>
                    <div style={{ width: '45%', marginLeft: '5%' }}>
                      <ScatterplotPca />
                    </div>
                  </div>
                  <h3>Another HTML element</h3>
                  <div style={{ width: '95%', height: '500px' }}>
                    <Heatmap />
                  </div>
                </div>
                <div style={{ width: '20%', height: '520px' }}>
                  <CellSets />
                  <GeneList />
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
    vw = vc.widget(page_esm=PAGE_ESM, page_mode=True, height=1100)
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
