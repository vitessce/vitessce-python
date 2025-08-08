import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Convert data manually
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        When running the Vitessce widget, data is converted on-the-fly into formats that the Vitessce JavaScript component can render.
        The converted data is stored in temporary files, preventing long-term use of the converted files.

        However, the data conversion utilities used by the widget are exposed so that their outputs can be saved to regular files.
        This allows the files to be saved for future use with the Vitessce web application by serving the files locally or moving the files onto an object storage system such as AWS S3 (for long-term storage and data sharing).

        This notebook demonstrates how to save the processed outputs of the `AnnDataWrapper` and `SnapWrapper` classes.
        """
    )
    return


@app.cell
def _():
    from vitessce import SnapWrapper
    from os.path import join
    from scipy.io import mmread
    import pandas as pd
    import numpy as np
    import json
    return SnapWrapper, join, json, mmread, pd


@app.cell
def _(join, mmread, pd):
    mtx = mmread(join('data', 'snapatac', 'filtered_cell_by_bin.mtx'))
    barcodes_df = pd.read_csv(join('data', 'snapatac', 'barcodes.txt'), header=None)
    bins_df = pd.read_csv(join('data', 'snapatac', 'bins.txt'), header=None)
    clusters_df = pd.read_csv(join('data', 'snapatac', 'umap_coords_clusters.csv'), index_col=0)
    return barcodes_df, bins_df, clusters_df, mtx


@app.cell
def _(SnapWrapper, barcodes_df, bins_df, clusters_df, join, json, mtx):
    zarr_filepath = join('data', 'snapatac', 'out.snap.multires.zarr')

    w = SnapWrapper(mtx, barcodes_df, bins_df, clusters_df)

    cells_json = w.create_cells_json()
    cell_sets_json = w.create_cell_sets_json()

    with open(join('data', 'snapatac', 'out.cells.json'), 'w') as f:
        json.dump(cells_json, f)

    with open(join('data', 'snapatac', 'out.cell-sets.json'), 'w') as f:
        json.dump(cell_sets_json, f)


    w.create_genomic_multivec_zarr(zarr_filepath)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
