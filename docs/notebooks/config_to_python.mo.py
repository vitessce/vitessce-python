import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Generate Python code to reconstruct a VitessceConfig instance""")
    return


@app.cell
def _():
    from vitessce import VitessceConfig, VitessceChainableConfig, VitessceConfigDatasetFile
    return (VitessceConfig,)


@app.cell
def _():
    from example_configs import dries as dries_config
    return (dries_config,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Load a view config from a dict representation""")
    return


@app.cell
def _(VitessceConfig, dries_config):
    vc = VitessceConfig.from_dict(dries_config)
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Print to JSON""")
    return


@app.cell
def _(vc):
    import json
    print(json.dumps(vc.to_dict(), indent=2))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Print to Python

    The `vc.to_python` function generates formatted Python code which can be used to re-generate the `vc` instance.
    """
    )
    return


@app.cell
def _(vc):
    imports, code = vc.to_python()
    return code, imports


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The first value returned is a list of classes used by the code snippet.""")
    return


@app.cell
def _(imports):
    imports
    return


@app.cell
def _(code):
    print(code)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The second value is the code snippet. When evaluated, the result will be a new `VitessceConfig` instance.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Evaluate the code and render a Vitessce widget""")
    return


@app.cell
def _(code):
    reconstructed_vc = eval(code)
    return (reconstructed_vc,)


@app.cell
def _(reconstructed_vc):
    reconstructed_vc.widget()
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
