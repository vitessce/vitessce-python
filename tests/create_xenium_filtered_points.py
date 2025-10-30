import os
from os.path import join, isfile, isdir
from urllib.request import urlretrieve
import zipfile
import shutil

# Used spatialdata==0.4.0 on October 30, 2025
from spatialdata import read_zarr, SpatialData


data_dir = "data"
zip_filepath = join(data_dir, "xenium_rep1_io.spatialdata.zarr.zip")
spatialdata_filepath = join(data_dir, "xenium_rep1_io.spatialdata.zarr")


if not isdir(spatialdata_filepath):
    if not isfile(zip_filepath):
        os.makedirs(data_dir, exist_ok=True)
        urlretrieve('https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep1_io.zip', zip_filepath)
    with zipfile.ZipFile(zip_filepath, "r") as zip_ref:
        zip_ref.extractall(data_dir)
        os.rename(join(data_dir, "data.zarr"), spatialdata_filepath)

        # This Xenium dataset has an AnnData "raw" element.
        # Reference: https://github.com/giovp/spatialdata-sandbox/issues/55
        raw_dir = join(spatialdata_filepath, "tables", "table", "raw")
        if isdir(raw_dir):
            shutil.rmtree(raw_dir)

sdata = read_zarr(spatialdata_filepath)

ddf = sdata.points["transcripts"]

# 2. Define a function to take every 100th row from a partition


def select_every_200th(partition):
    # Each 'partition' is a Pandas DataFrame
    # .iloc[::100] is the efficient pandas way to get every 100th row
    return partition.iloc[::200]


# 3. Apply this function to every partition in the Dask DataFrame
result = ddf.map_partitions(select_every_200th)

# 4. Compute the result to see it
filtered_ddf = result[["x", "y", "z", "feature_name", "cell_id"]]

small_sdata = SpatialData(points={"transcripts": filtered_ddf})

small_sdata.write("xenium_rep1_io.points_only.spatialdata.zarr", overwrite=True)
