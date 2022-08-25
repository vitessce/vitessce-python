import numpy as np
import zarr
from ome_zarr.writer import write_image

def rgb_img_to_ome_zarr(img_arr, output_path, img_name="Image", chunks=(1, 256, 256), axes="cyx"):
    img_arr = img_arr.astype(np.dtype('uint8'))

    default_window = {
        "start": 0,
        "min": 0,
        "max": 255,
        "end": 255
    }

    z_root = zarr.open_group(output_path)
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