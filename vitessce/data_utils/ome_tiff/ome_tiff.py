import numpy as np
from tifffile import TiffWriter


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

    tiff_writer = TiffWriter(output_path, ome=True)
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
    tiff_writer = TiffWriter(output_path, ome=True)
    tiff_writer.write(
        img_arr,
        metadata={
            'axes': axes,
            'Channel': {'Name': channel_names},
        }
    )
    tiff_writer.close()
