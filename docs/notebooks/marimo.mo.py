import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from vitessce import VitessceConfig
    from example_configs import dries as dries_config
    return VitessceConfig, dries_config


@app.cell
def _(VitessceConfig, dries_config, mo):
    vc = VitessceConfig.from_dict(dries_config)
    vw = mo.ui.anywidget(vc.widget())
    vw
    return (vw,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""As you interact with the widget above, watch as the Vitessce configuration below dynamically updates. For instance, as you zoom in the scatterplot, check the value of `coordinationSpace.embeddingZoom`""")
    return


@app.cell
def _(vw):
    vw._config
    return


@app.cell
def _(vw):
    vw._config["coordinationSpace"]["embeddingZoom"]
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
