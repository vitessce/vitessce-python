import numpy as np
import zarr
from ome_zarr.writer import write_image
from .anndata import cast_arr


def rgb_img_to_ome_zarr(img_arr, output_path, img_name="Image", chunks=(1, 256, 256), axes="cyx"):
    img_arr = img_arr.astype(np.dtype('uint8'))

    default_window = {
        "start": 0,
        "min": 0,
        "max": 255,
        "end": 255
    }

    z_root = zarr.open_group(output_path, mode="w")

    write_image(
        image=img_arr,
        group=z_root,
        axes=axes,
        omero={
            "name": img_name,
            "version": "0.3",
            "rdefs": {},
            "channels": [
                {
                    "label": "R",
                    "color": "FF0000",
                    "window": default_window
                },
                {
                    "label": "G",
                    "color": "00FF00",
                    "window": default_window
                },
                {
                    "label": "B",
                    "color": "0000FF",
                    "window": default_window
                }
            ]
        },
        chunks=chunks
    )


def multiplex_img_to_ome_zarr(img_arr, channel_names, output_path, img_name="Image", chunks=(1, 256, 256), axes="cyx", channel_colors=None):
    img_arr = cast_arr(img_arr)

    dtype_info = np.iinfo(img_arr.dtype) if img_arr.dtype.kind == 'u' or img_arr.dtype.kind == 'i' else np.finfo(img_arr.dtype)

    default_window = {
        "start": 0,
        "min": 0,
        "max": dtype_info.max,
        "end": dtype_info.max
    }

    z_root = zarr.open_group(output_path, mode="w")

    write_image(
        image=img_arr,
        group=z_root,
        axes=axes,
        omero={
            "name": img_name,
            "version": "0.3",
            "rdefs": {},
            "channels": [
                {
                    "label": channel_name,
                    "color": channel_colors[channel_name] if channel_colors is not None else "FFFFFF",
                    "window": default_window
                }
                for channel_name
                in channel_names
            ]
        },
        chunks=chunks
    )
