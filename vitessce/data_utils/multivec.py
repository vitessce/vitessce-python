import math
import zarr
import numpy as np
import pandas as pd
from tqdm import tqdm

from .anndata import to_dense
from .entities import GenomicProfiles


def adata_to_multivec_zarr(adata, output_path, obs_set_col, obs_set_name, obs_set_vals=None, var_interval_col="interval", layer_key=None, assembly="hg38", starting_resolution=5000, chr_subset=None):
    """
    Convert an AnnData object containing a cell-by-bin matrix to a Multivec-Zarr store.

    :param adata: The object to convert.
    :type adata: anndata.AnnData
    :param output_path: The path to the output Zarr store.
    :type output_path: str
    :param obs_set_col: The name of the column in adata.obs that contains the cluster IDs.
    :type obs_set_col: str
    :param obs_set_name: The name of the cluster set.
    :type obs_set_name: str
    :param obs_set_vals: The cluster IDs to include in the output. If None, all cluster IDs will be included. This can be used to override the order of the cluster IDs.
    :type obs_set_vals: list[str] or None
    :param var_interval_col: The name of the column in adata.var that contains the bin interval strings. By default, "interval".
    :type var_interval_col: str
    :param layer_key: The name of the layer in adata.layers to use. If None, adata.X will be used. By default, None.
    :type layer_key: str or None
    :param assembly: The name of the genome assembly. By default, "hg38".
    :type assembly: str
    :param starting_resolution: The starting resolution of the data. By default, 5000.
    :type starting_resolution: int
    :param chr_subset: For debugging purposes, a subset of chromosomes to process. If None, all chromosomes in the assembly will be processed. By default, None.
    :type chr_subset: list[str] or None
    """
    in_mtx = adata.layers[layer_key] if layer_key is not None else adata.X
    in_barcodes_df = adata.obs
    in_bins_df = adata.var

    # Ensure that in_bins_df has a sequential integer index
    in_bins_df = in_bins_df.reset_index()

    in_mtx = to_dense(in_mtx)  # TODO: is this necessary?

    # The bin datafram consists of one column like chrName:binStart-binEnd
    def convert_bin_name_to_chr_name(bin_name):
        try:
            return bin_name[:bin_name.index(':')]
        except ValueError:
            return np.nan

    def convert_bin_name_to_chr_start(bin_name):
        try:
            return int(bin_name[bin_name.index(':') + 1:bin_name.index('-')])
        except ValueError:
            return np.nan

    def convert_bin_name_to_chr_end(bin_name):
        try:
            return int(bin_name[bin_name.index('-') + 1:])
        except ValueError:
            return np.nan

    # Keep only the interval column
    in_bins_df = in_bins_df[[var_interval_col]]
    in_bins_df = in_bins_df.rename(columns={var_interval_col: "interval"})
    in_bins_df["chr_name"] = in_bins_df["interval"].apply(
        convert_bin_name_to_chr_name)
    in_bins_df["chr_start"] = in_bins_df["interval"].apply(
        convert_bin_name_to_chr_start)
    in_bins_df["chr_end"] = in_bins_df["interval"].apply(
        convert_bin_name_to_chr_end)

    # Drop any rows that had incorrect bin strings (missing a chromosome name, bin start, or bin end value).
    in_bins_df = in_bins_df.dropna(
        subset=["chr_name", "chr_start", "chr_end"]).copy()

    # Ensure that the columns have the expected types.
    in_bins_df["chr_name"] = in_bins_df["chr_name"].astype(str)
    in_bins_df["chr_start"] = in_bins_df["chr_start"].astype(int)
    in_bins_df["chr_end"] = in_bins_df["chr_end"].astype(int)

    in_bins_df = in_bins_df.reset_index(drop=True)

    interval_sizes = in_bins_df.apply(lambda row: row["chr_end"] - row["chr_start"], axis='columns')
    max_interval = interval_sizes.max()
    if max_interval > starting_resolution:
        raise ValueError("Starting resolution is smaller than largest interval.")

    # Round bins
    in_bins_df["chr_start_round"] = in_bins_df["chr_start"].apply(lambda x: math.floor(x / starting_resolution) * starting_resolution + 1)
    in_bins_df["chr_end_round"] = in_bins_df["chr_start_round"].apply(lambda x: x + starting_resolution - 1)
    # TODO: do the values need to be scaled based on the ratio of the original size of the interval to the rounded size?

    # Replace the original start/end values
    in_bins_df["chr_start"] = in_bins_df["chr_start_round"]
    in_bins_df["chr_end"] = in_bins_df["chr_end_round"]
    in_bins_df = in_bins_df.drop(columns=["chr_start_round", "chr_end_round"])
    in_bins_df["interval"] = in_bins_df.apply(lambda r: f"{r['chr_name']}:{r['chr_start']}-{r['chr_end']}", axis='columns')

    # Use provided obs_set_vals if present, since these may be ordered
    # in a particular way.
    cluster_ids = in_barcodes_df[obs_set_col].unique().tolist() if obs_set_vals is None else obs_set_vals

    cluster_paths = [[obs_set_name, cluster_id] for cluster_id in cluster_ids]

    # Create the Zarr store for the outputs.
    out_f = zarr.open(output_path, mode='w')

    genomic_profiles = GenomicProfiles(
        out_f, profile_paths=cluster_paths, assembly=assembly, starting_resolution=starting_resolution
    )
    chrom_name_to_length = genomic_profiles.chrom_name_to_length

    # For debugging purposes, allow the user to specify a subset of chromosomes to process.
    chrom_names = chr_subset if chr_subset is not None else list(chrom_name_to_length.keys())

    # Create each chromosome dataset.
    for chr_name in tqdm(chrom_names):
        chr_len = chrom_name_to_length[chr_name]
        # The bins dataframe frustratingly does not contain every bin.
        # We need to figure out which bins are missing.

        # We want to check for missing bins in each chromosome separately,
        # otherwise too much memory is used during the join step.
        chr_bins_in_df = in_bins_df.loc[in_bins_df["chr_name"] == chr_name]
        if chr_bins_in_df.shape[0] == 0:
            print("Warning: No bins found for chromosome", chr_name)
            # No processing or output is necessary if there is no data for this chromosome.
            # Continue on through all resolutions of this chromosome to the next chromosome.
            continue
        # Determine the indices of the matrix at which the bins for this chromosome start and end.
        chr_bin_i_start = int(chr_bins_in_df.head(1).iloc[0].name)
        # +1 because the end index is exclusive.
        chr_bin_i_end = int(chr_bins_in_df.tail(1).iloc[0].name) + 1

        # Extract the part of the matrix corresponding to the current chromosome.
        chr_mtx = in_mtx[:, chr_bin_i_start:chr_bin_i_end]

        # Create a list of the "ground truth" bins (all bins from position 0 to the end of the chromosome).
        # We will join the input bins onto this dataframe to determine which bins are missing.
        chr_bins_gt_df = pd.DataFrame()
        chr_bins_gt_df["chr_start"] = np.arange(0, math.ceil(
            chr_len / starting_resolution)) * starting_resolution
        chr_bins_gt_df["chr_end"] = chr_bins_gt_df["chr_start"] + \
            starting_resolution
        chr_bins_gt_df["chr_start"] = chr_bins_gt_df["chr_start"] + 1
        chr_bins_gt_df["chr_start"] = chr_bins_gt_df["chr_start"].astype(
            int)
        chr_bins_gt_df["chr_end"] = chr_bins_gt_df["chr_end"].astype(int)
        chr_bins_gt_df["chr_name"] = chr_name
        chr_bins_gt_df[0] = chr_bins_gt_df.apply(lambda r: f"{r['chr_name']}:{r['chr_start']}-{r['chr_end']}", axis='columns')

        # We will add a new column "i", which should match the _old_ index, so that we will be able join with the data matrix on the original indices.
        # For the new rows, we will add values for the "i" column that are greater than any of the original indices,
        # to prevent any joining with the incoming data matrix onto these bins for which the data is missing.
        chr_bins_in_df = chr_bins_in_df.reset_index(drop=True)
        chr_bins_in_df["i"] = chr_bins_in_df.index.values
        chr_bins_gt_df["i"] = chr_bins_gt_df.index.values + \
            (in_mtx.shape[1] + 1)

        # Set the full bin string column as the index of both data frames.
        chr_bins_gt_df = chr_bins_gt_df.set_index(0)
        chr_bins_in_df = chr_bins_in_df.set_index("interval")

        # Join the input bin subset dataframe right onto the full bin ground truth dataframe.
        chr_bins_in_join_df = chr_bins_in_df.join(
            chr_bins_gt_df, how='right', lsuffix="", rsuffix="_gt")
        # The bins which were not present in the input will have NaN values in the "i" column.
        # For these rows, we replace the NaN values with the much higher "i_gt" values which will not match to any index of the data matrix.
        chr_bins_in_join_df["i"] = chr_bins_in_join_df.apply(
            lambda r: r['i'] if pd.notna(r['i']) else r['i_gt'], axis='columns').astype(int)

        # Clean up the joined data frame by removing unnecessary columns.
        chr_bins_in_join_df = chr_bins_in_join_df.drop(
            columns=['chr_name', 'chr_start', 'chr_end', 'i_gt'])
        chr_bins_in_join_df = chr_bins_in_join_df.rename(
            columns={'chr_name_gt': 'chr_name', 'chr_start_gt': 'chr_start', 'chr_end_gt': 'chr_end'})

        # Create a dataframe from the data matrix, so that we can join to the joined bins dataframe.
        chr_mtx_df = pd.DataFrame(data=chr_mtx.T)

        chr_bins_i_df = chr_bins_in_join_df.drop(
            columns=['chr_name', 'chr_start', 'chr_end'])

        # Join the data matrix dataframe and the bins dataframe.
        # Bins that are missing from the data matrix will have "i" values higher than any of the data matrix dataframe row indices,
        # and therefore the data values for these bins in the resulting joined dataframe will all be NaN.
        chr_mtx_join_df = chr_bins_i_df.join(
            chr_mtx_df, how='left', on='i')
        # We fill in these NaN values with 0.
        chr_mtx_join_df = chr_mtx_join_df.fillna(value=0.0)

        # Drop the "i" column, since it is not necessary now that we have done the join.
        chr_mtx_join_df = chr_mtx_join_df.drop(columns=['i'])
        # Obtain the new full data matrix, which contains values for all bins of the chromosome.
        chr_mtx = chr_mtx_join_df.values.T

        # Fill in the Zarr store with data for each cluster.
        for cluster_index, cluster_id in enumerate(cluster_ids):
            # Get the list of cells in the current cluster.
            cluster_df = in_barcodes_df.loc[in_barcodes_df[obs_set_col]
                                            == cluster_id]
            cluster_cell_ids = cluster_df.index.values.tolist()
            cluster_cells_tf = (
                in_barcodes_df.index.to_series().isin(cluster_cell_ids)).values

            # Get the rows of the data matrix corresponding to the cells in this cluster.
            cluster_cell_by_bin_mtx = chr_mtx[cluster_cells_tf, :]
            # Take the sum of this cluster along the cells axis.
            cluster_profile = cluster_cell_by_bin_mtx.sum(axis=0)

            # For some reason the matrix can contain intervals past the end of the
            # chromosome according to the length from negspy,
            # so we only keep those bins that fit.
            profile_len = math.ceil(chr_len / starting_resolution)
            # TODO: raise warning if the cluster_profile length is longer than profile_len?

            genomic_profiles.add_profile(
                cluster_profile[0:profile_len], chr_name, cluster_index)
