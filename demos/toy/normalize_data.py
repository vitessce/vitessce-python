import spatialdata as sd
import scanpy as sc

# Read in the data as a SpatialData object
sdata = sd.read_zarr(snakemake.input[0])

# Get the AnnData object that the SpatialData object contains
adata = sdata.table

# The ScanPy package is compatible with AnnData objects
# References:
# - https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# - https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)

sdata.table = adata

# Write the normalized data back to the SpatialData object
sdata.write_zarr(snakemake.output[0])
