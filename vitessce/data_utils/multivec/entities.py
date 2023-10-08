import negspy.coordinates as nc
import numpy as np
import math


class ArgumentLengthDoesNotMatchCellIdsException(Exception):
    pass


class NodeNotFoundException(Exception):
    pass


class GenomicProfiles():

    """
    Generic class for representing genomic profiles.
    """

    def __init__(self, f, profile_paths, assembly='hg38', starting_resolution=5000, name="Genomic Profiles"):
        """
        Constructor method

        :param f: The opened Zarr store object.
        :type f: zarr.Group
        :param list[list[str]] profile_paths: A list of cell set paths, one path for each profile.
        :param str assembly: The genome assembly to use for chromosome lengths, passed to negspy. By default, 'hg38'.
        :param int starting_resolution: The starting resolution. By default, 5000.
        :param str name: The name for this set of profiles. By default, 'Genomic Profiles'.
        """

        self.f = f

        num_profiles = len(profile_paths)

        compressor = 'default'

        chromosomes = [str(chr_name) for chr_name in nc.get_chromorder(
            assembly)[:25]]  # TODO: should more than chr1-chrM be used?
        chroms_length_arr = np.array(
            [nc.get_chrominfo(assembly).chrom_lengths[x] for x in chromosomes], dtype="i8")
        chroms_cumsum_arr = np.concatenate(
            (np.array([0]), np.cumsum(chroms_length_arr)))

        chrom_name_to_length = dict(zip(chromosomes, chroms_length_arr))
        chrom_name_to_cumsum = dict(zip(chromosomes, chroms_cumsum_arr))

        # Prepare to fill in resolutions datasets.
        resolutions = [starting_resolution * (2 ** x) for x in range(16)]

        chromosomes_group = f.create_group("chromosomes")
        for chr_name, chr_len in chrom_name_to_length.items():
            chr_group = chromosomes_group.create_group(chr_name)
            # Create each resolution group.
            for resolution in resolutions:
                chr_shape = (num_profiles, math.ceil(chr_len / resolution))
                chr_group.create_dataset(str(
                    resolution), shape=chr_shape, dtype="f4", fill_value=np.nan, compressor=compressor)

        # f.attrs should contain the properties required for HiGlass's "tileset_info" requests.
        f.attrs['row_infos'] = [
            {"path": profile_path}
            for profile_path in profile_paths
        ]
        f.attrs['resolutions'] = sorted(resolutions, reverse=True)
        f.attrs['shape'] = [num_profiles, 256]
        f.attrs['name'] = name
        f.attrs['coordSystem'] = assembly

        self.resolutions = resolutions
        self.chromosomes = chromosomes
        self.chromosomes_group = chromosomes_group
        self.chrom_name_to_length = chrom_name_to_length
        self.num_profiles = num_profiles

        # https://github.com/zarr-developers/zarr-specs/issues/50
        f.attrs['multiscales'] = [
            {
                "version": "0.1",
                "name": chr_name,
                "datasets": [
                    {"path": f"chromosomes/{chr_name}/{resolution}"}
                    for resolution in sorted(resolutions, reverse=True)
                ],
                "type": "zarr-multivec",
                "metadata": {
                    "chromoffset": int(chrom_name_to_cumsum[chr_name]),
                    "chromsize": int(chr_len),
                }
            }
            for (chr_name, chr_len) in list(zip(chromosomes, chroms_length_arr))
        ]

    def add_profile(self, values, chr_name, profile_index):
        """
        Add a single genomic profile to the output store. This function will aggregate for each resolution.

        :param values: A profile array for one chromosome.
        :type values: np.array
        :param str chr_name: The name the chromosome corresponding to this array.
        :param int profile_index: The index of this profile among the list of profiles.
        """
        chromosomes_group = self.chromosomes_group
        resolutions = self.resolutions
        resolution_exps = [(2**x) for x in range(len(resolutions))]

        chr_len = self.chrom_name_to_length[chr_name]
        # Fill in the data for this cluster and chromosome at each resolution.
        for resolution, resolution_exp in zip(resolutions, resolution_exps):
            arr_len = math.ceil(chr_len / resolution)

            # Pad the array of values with zeros if necessary before reshaping.
            padding_len = resolution_exp - (values.shape[0] % resolution_exp)
            if values.shape[0] % resolution_exp > 0:
                values = np.concatenate((values, np.zeros((padding_len,))))
            # Reshape to be able to sum every `resolution_exp` number of values.
            arr = np.reshape(values, (-1, resolution_exp)).sum(axis=-1)

            padding_len = arr_len - arr.shape[0]
            if padding_len > 0:
                arr = np.concatenate((arr, np.zeros((padding_len,))))
            # Set the array in the Zarr store.
            chromosomes_group[chr_name][str(
                resolution)][profile_index, :] = arr
