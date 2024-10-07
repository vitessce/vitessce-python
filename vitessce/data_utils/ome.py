import numpy as np
import zarr
from ome_zarr.writer import write_image
from tifffile import TiffWriter
from .anndata import cast_arr


def needs_bigtiff(img_arr_shape):
    """
    Helper function to determine if an image array is too large for standard TIFF format.

    :param img_arr_shape: The shape of the image array.
    :type img_arr_shape: tuple[int]
    :return: True if the image array is too large for standard TIFF format, False otherwise.
    :rtype: bool
    """
    num_pixels = 1
    for n in img_arr_shape.shape:
        num_pixels *= n
    return (num_pixels > 2**32)


def rgb_img_to_ome_tiff(img_arr, output_path, img_name="Image", axes="CYX"):
    """
    Convert an RGB image to OME-TIFF.

    :param img_arr: The image as a 3D array.
    :type img_arr: np.array
    :param str output_path: The path to save the Zarr store.
    :param str img_name: The name of the image to include in the omero.name NGFF metadata field.
    :param str axes: The array axis ordering. By default, "CYX"
    """
    img_arr = img_arr.astype(np.dtype('uint8'))
    bigtiff = needs_bigtiff(img_arr.shape)

    tiff_writer = TiffWriter(output_path, ome=True, bigtiff=bigtiff)
    tiff_writer.write(
        img_arr,
        metadata={
            'axes': axes,
            'Channel': {'Name': ['R', 'G', 'B']},
        }
    )
    tiff_writer.close()


def multiplex_img_to_ome_tiff(img_arr, channel_names, output_path, axes="CYX"):
    """
    Convert a multiplexed image to OME-TIFF.

    :param img_arr: The image as a 3D, 4D, or 5D array.
    :type img_arr: np.array
    :param list[str] channel_names: A list of channel names to include in the omero.channels[].label NGFF metadata field.
    :param str output_path: The path to save the Zarr store.
    :param str axes: The array axis ordering. By default, "CYX"
    """
    bigtiff = needs_bigtiff(img_arr.shape)

    tiff_writer = TiffWriter(output_path, ome=True, bigtiff=bigtiff)
    tiff_writer.write(
        img_arr,
        metadata={
            'axes': axes,
            'Channel': {'Name': channel_names},
        }
    )
    tiff_writer.close()


def rgb_img_to_ome_zarr(img_arr, output_path, img_name="Image", chunks=(1, 256, 256), axes="cyx", **kwargs):
    """
    Convert an RGB image to OME-Zarr v0.3.

    :param img_arr: The image as a 3D array.
    :type img_arr: np.array
    :param str output_path: The path to save the Zarr store.
    :param str img_name: The name of the image to include in the omero.name NGFF metadata field.
    :param tuple[int] chunks: The chunk sizes of each axis. By default, (1, 256, 256).
    :param str axes: The array axis ordering. By default, "cyx"
    """
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
        storage_options=dict(chunks=chunks),
        **kwargs,
    )
    z_root.attrs["omero"] = {
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
    }


def multiplex_img_to_ome_zarr(img_arr, channel_names, output_path, img_name="Image", chunks=(1, 256, 256), axes="cyx", channel_colors=None):
    """
    Convert a multiplexed image to OME-Zarr v0.3.

    :param img_arr: The image as a 3D, 4D, or 5D array.
    :type img_arr: np.array
    :param list[str] channel_names: A list of channel names to include in the omero.channels[].label NGFF metadata field.
    :param str output_path: The path to save the Zarr store.
    :param str img_name: The name of the image to include in the omero.name NGFF metadata field.
    :param tuple[int] chunks: The chunk sizes of each axis. By default, (1, 256, 256).
    :param str axes: The array axis ordering. By default, "cyx"
    :param channel_colors: Dict mapping channel names to color strings to use for the omero.channels[].color NGFF metadata field. If provided, keys should match channel_names. By default, None to use "FFFFFF" for all channels.
    :type channel_colors: dict or None
    """
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
        storage_options=dict(chunks=chunks)
    )
    z_root.attrs["omero"] = {
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
    }
