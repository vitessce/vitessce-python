import negspy.coordinates as nc
import numpy as np
import math


class ArgumentLengthDoesNotMatchCellIdsException(Exception):
    pass


class NodeNotFoundException(Exception):
    pass


class Cells:

    """
    Generic Cells class for constructing the json needed for client side rendering of cell segmentations/scatterplots (UMAP, PCA etc.).

    :param json The json resulting from various calls to add_mapping, add_polygon_outline etc. that can be served to the client.
    """

    def __init__(self, cell_ids=[]):
        """
        Constructor method

        :param list cell_ids: A list of cell ids to be shown in Vitessce.  The order of these will be used to determine the order of future additions to this class, like segmentations and scatterplot coordinates.
        """
        self._cell_ids = cell_ids
        self.json = dict(zip(cell_ids, [{} for _ in cell_ids]))

    def add_mapping(self, name, coords):
        """
        Add a (dimensionality reduction) scatterplot mapping to each cell.

        :param str name: The unique identifier for the mapping, like UMAP, tSNE or PCA.
        :param list coords: A list of lists like [[1, 2], [3, 4], ...] in the order of cell_ids for each cell to be mapped to a scatterplot coorindate.
        """
        if len(coords) != len(self._cell_ids):
            raise ArgumentLengthDoesNotMatchCellIdsException(
                'Coordinates length does not match Cell IDs Length')
        if not isinstance(name, str):
            raise TypeError(
                'name argument needs to be a string for adding a scatterplot mapping')
        for idx, id in enumerate(self._cell_ids):
            if 'mappings' not in self.json[id]:
                self.json[id]['mappings'] = {name: coords[idx]}
            else:
                self.json[id]['mappings'][name] = coords[idx]

    def add_centroids(self, centroids):
        """
        Add a centroid for a spatial segmentation outline to each cell.

        :param list centroids: A list of lists like [[1, 2], [3, 4], ...] in the order of cell_ids for each cell to be mapped to a centroid coorindate.
        """
        if len(centroids) != len(self._cell_ids):
            raise ArgumentLengthDoesNotMatchCellIdsException(
                'Centroid length does not match Cell IDs Length')
        if not isinstance(centroids, list) or any([len(centroid) != 2 or not isinstance(centroid, list) for centroid in centroids]):
            raise TypeError('Centroids should be a list of two element lists')
        for idx, id in enumerate(self._cell_ids):
            self.json[id]['xy'] = centroids[idx]

    def add_polygon_outline(self, polygon_outline):
        """
        Add a polygon for a spatial segmentation outline to each cell.

        :param list polygon_outline: A list of lists of lists like [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]...] in the order of cell_ids for each cell to be mapped to its segmentation.
        """
        if len(polygon_outline) != len(self._cell_ids):
            raise ArgumentLengthDoesNotMatchCellIdsException(
                'Segmentations length does not match Cell IDs Length')
        for idx, id in enumerate(self._cell_ids):
            if not isinstance(polygon_outline[idx], list) or any([len(coord) != 2 or not isinstance(coord, list) for coord in polygon_outline[idx]]):
                raise TypeError(
                    f'Polygon outline for {id} should be a list of two element lists i.e xy coordinates')
            self.json[id]['poly'] = polygon_outline[idx]


class CellSets:

    """
    Generic CellSets class for constructing the json needed for client side rendering of the cell sets.

    :param json The json resulting from various calls to add_node that can be served to the client.
    """

    def __init__(self):
        """
        Constructor method
        """

        self.json = {
            "datatype": "cell",
            "version": "0.1.2",
            "tree": []
        }

    def add_level_zero_node(self, name):
        """
        Add a new level zero node to the root of the tree.

        :param str name: Name for the new node
        """
        self.json['tree'].append({
            "name": name,
            "children": []
        })

    def add_node(self, name, parent_path, cell_set=None):
        """
        Add a node to a parent node.

        :param str name: Name for the new node
        :param list parent_path: List of strings representing the internal nodes to traverse to reach the desired parent node to which we will add the new node, like ['epithelial', 'meso-epithelial']
        :param list cell_set: List of cell ids which will be added to the new node as part of the set.
        """
        parent_node = self._tree_find_node_by_path(parent_path)
        if parent_node is None:
            raise NodeNotFoundException(
                f'No node with path {parent_path} found to add {name} to')
        new_node = {"name": name}
        if cell_set:
            new_node['set'] = cell_set
        if 'children' not in parent_node:
            parent_node['children'] = [new_node]
        else:
            parent_node['children'].append(new_node)

    def _find_node_by_path(self, node, path, curr_index):
        curr_node_name = path[curr_index]
        if node['name'] == curr_node_name:
            if curr_index == len(path) - 1:
                return node
            if 'children' in node:
                found_nodes = [
                    self._find_node_by_path(child, path, curr_index + 1) for child in node['children']
                ]
                found_nodes_not_none = [n for n in found_nodes if n]
                if len(found_nodes_not_none) == 1:
                    return found_nodes[0]
        return None

    def _tree_find_node_by_path(self, path):
        found_nodes = [self._find_node_by_path(
            node, path, 0) for node in self.json['tree']]
        found_nodes_not_none = [n for n in found_nodes if n]
        if len(found_nodes_not_none) == 1:
            return found_nodes_not_none[0]
        return None


class Molecules():

    """
    Generic Molecules class for constructing the json needed for client side rendering of spot data.

    :param json The json resulting from various calls to add_molecule.
    """

    def __init__(self):
        """
        Constructor method
        """
        self.json = {}

    def add_molecule(self, name, coords):
        """
        Add a molecules to a parent node.

        :param str name: Name for the new molecules
        :param list coords: A list of lists like [[1, 2], [3, 4], ...] or [[1, 2, 3], [3, 4, 5], ...] which denote where in xy space the spot data should be placed for the desired name.
        """
        self.json[name] = coords


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
