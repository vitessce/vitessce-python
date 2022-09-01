import argparse
import h5py
import numpy as np
from vitessce.data_utils import (
    multiplex_img_to_ome_zarr,
)


def convert_to_ome_zarr(args):
    with h5py.File(args.input, "r") as f:
        channel_names = list(f.keys())
        first_channel = channel_names[0]
        img_shape = (len(channel_names), *f[first_channel].shape)
        img_dtype = f[first_channel].dtype
        img_arr = np.zeros(img_shape, dtype=img_dtype)
        for channel_i, channel_name in enumerate(channel_names):
            img_arr[channel_i, :, :] = f[channel_name]

        img_arr = img_arr.astype(np.dtype("<u2"))
        img_arr = np.transpose(img_arr, axes=(0, 2, 1))

        multiplex_img_to_ome_zarr(
            img_arr,
            channel_names,
            args.output,
            img_name="Image",
            chunks=(1, 512, 512),
            axes="cyx",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input HDF5 file"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output OME Zarr store"
    )
    args = parser.parse_args()
    convert_to_ome_zarr(args)
