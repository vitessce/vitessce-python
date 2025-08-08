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
        # Using an existing view config dict
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 1. Import dependencies

        We need to import the `VitessceConfig` class.
        """
    )
    return


@app.cell
def _():
    from vitessce import VitessceConfig
    return (VitessceConfig,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Import config of interest as a dict

        The Vitessce config at its lowest level is a JSON object or a Python `dict` that specifies a layout for the Vitessce components will be rendered in the widget. These components may be scatterplots, spatial plots, heatmaps, or control components. The config also describes the datasets and files that will be visualized.

        The `vitessce` package provides helper functions and classes to simplify the process of defining Vitessce configs. Those functions are demonstrated in the other notebooks. The helper functions are intended to make visualization of local datasets easy. However, in this case, we are importing a config that has been pre-defined in the file `example_configs.py`, in which the dataset being visualized is stored remotely in AWS S3 (rather than locally).
        """
    )
    return


@app.cell
def _():
    from example_configs import dries as dries_config
    return (dries_config,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 3. Create the Vitessce widget

        Create the widget by creating a new config instance using the `VitessceConfig.from_dict` static method, then calling the config's `.widget()` method.

        Render the widget by placing the widget variable on its own line at the end of the notebook cell.
        """
    )
    return


@app.cell
def _(VitessceConfig, dries_config):
    vc = VitessceConfig.from_dict(dries_config)
    vw = vc.widget()
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
