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
        # Visualization of AnnData object containing an image in `uns`
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Note: This approach to storing images within AnnData objects is no longer recommended now that [SpatialData](https://spatialdata.scverse.org/en/stable/) has been introduced.
        """
    )
    return


@app.cell
def _():
    import scanpy as sc
    import numpy as np
    from vitessce.data_utils import rgb_img_to_ome_zarr, VAR_CHUNK_SIZE
    from vitessce import (
        VitessceConfig,
        AnnDataWrapper,
        ImageOmeZarrWrapper,
    )
    from os.path import join
    return (
        AnnDataWrapper,
        ImageOmeZarrWrapper,
        VAR_CHUNK_SIZE,
        VitessceConfig,
        join,
        np,
        rgb_img_to_ome_zarr,
        sc,
    )


@app.cell
def _(join):
    output_img = join("data", "V1_Human_Lymph_Node.ome.zarr")
    output_adata = join("data", "V1_Human_Lymph_Node.anndata.zarr")
    return output_adata, output_img


@app.cell
def _(VAR_CHUNK_SIZE, np, output_adata, output_img, rgb_img_to_ome_zarr, sc):
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node", include_hires_tiff=True)

    # Write img_arr to OME-Zarr.
    # Need to convert images from interleaved to non-interleaved (color axis should be first).
    img_hires = adata.uns['spatial']['V1_Human_Lymph_Node']['images']['hires']
    img_arr = np.transpose(img_hires, (2, 0, 1))
    # Convert values from [0, 1] to [0, 255].
    img_arr *= 255.0

    # First, save the image to an OME-Zarr image format
    rgb_img_to_ome_zarr(img_arr, output_img, axes="cyx", chunks=(1, 256, 256), img_name="Image")
    # Second, save the AnnData object to Zarr format
    adata.write_zarr(output_adata, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
    return


@app.cell
def _(
    AnnDataWrapper,
    ImageOmeZarrWrapper,
    VitessceConfig,
    output_adata,
    output_img,
):
    vc = VitessceConfig(schema_version="1.0.17", name="AnnData with image")
    dataset = vc.add_dataset("My dataset").add_object(
        AnnDataWrapper(
            adata_path=output_adata,
        
        )
    ).add_object(
        ImageOmeZarrWrapper(
            img_path=output_img,
        )
    )

    spatial_view = vc.add_view("spatialBeta", dataset=dataset)
    lc_view = vc.add_view("layerControllerBeta", dataset=dataset)

    vc.layout(spatial_view | lc_view);
    return (vc,)


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
