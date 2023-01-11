from constants import COLUMNS, CL_ROOT_ID
from utils import (
    load_cl_obo_graph,
    init_cell_sets_tree,
    dict_to_tree,
    sort_paths_up_cell_ontology,
    get_paths_up_cell_ontology,
    fill_in_dict_from_path,
    remove_any_from_dict_levels,
)


def generate_leiden_cluster_cell_sets(df):
    """
    Generate a tree of cell sets representing
    the clusters from the Leiden clustering
    algorithm.
    """
    tree = init_cell_sets_tree()

    leiden_clusters_children = []
    for cluster_name, cluster_df in df.groupby("leiden"):
        leiden_clusters_children.append({
            "name": cluster_name,
            "set": [
                [x, None]
                for x in cluster_df[COLUMNS.CELL_ID.value].unique().tolist()
            ],
        })

    tree["tree"].append({
        "name": "Leiden Clustering",
        "children": leiden_clusters_children
    })

    return tree


def generate_cell_type_flat_cell_sets(df):
    """
    Generate a tree of cell sets
    representing the cell type annotations,
    arranged on one level (not heirarchical).
    """
    tree = init_cell_sets_tree()

    cell_type_annotation_children = []
    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):
        set_cell_ids = cell_type_df[COLUMNS.CELL_ID.value].values.tolist()
        set_cell_scores = (
            cell_type_df[COLUMNS.PREDICTION_SCORE.value]
            .values.tolist()
        )
        set_value = [list(x) for x in zip(set_cell_ids, set_cell_scores)]

        cell_type_annotation_children.append({
            "name": cell_type,
            "set": set_value,
        })

    tree["tree"].append({
        "name": "Cell Type Annotations",
        "children": cell_type_annotation_children
    })
    return tree


def generate_cell_type_cell_sets(df, cl_obo_file):
    """
    Generate a tree of cell sets
    for hierarchical cell type annotations.
    """
    tree = init_cell_sets_tree()

    # Load the cell ontology DAG
    graph, id_to_name, name_to_id = load_cl_obo_graph(cl_obo_file)

    ancestors_and_sets = []

    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):
        try:
            node_id = name_to_id[cell_type]
        except KeyError:
            print(
                f"ERROR: annotation '{cell_type}' does "
                "not match any node in the cell ontology."
            )
            continue

        # Get all of the possible paths up to the root
        # from the current node.
        paths_up = get_paths_up_cell_ontology(graph, node_id)
        # Get the names of each node in each path.
        named_paths_up = [
            [id_to_name[n_id] for n_id in path_nodes]
            for path_nodes in paths_up
        ]
        print(
            f"WARNING: {id_to_name[node_id]} has {len(paths_up)} paths"
            f" up to {CL_ROOT_ID} ({id_to_name[CL_ROOT_ID]})."
        )

        # Sort potential paths "up the hierarchy" by our preferences,
        # to avoid "functional" parent nodes like "motile cell"
        sorted_named_paths_up = sort_paths_up_cell_ontology(named_paths_up)

        named_ancestors = sorted_named_paths_up[0]
        named_ancestors_reversed = list(reversed(named_ancestors))
        # Get a list of (cell_id, prediction_score) tuples for the set.
        set_value = [
            list(x)
            for x in zip(
                cell_type_df[COLUMNS.CELL_ID.value].values.tolist(),
                cell_type_df[COLUMNS.PREDICTION_SCORE.value].values.tolist()
            )
        ]

        ancestors_and_sets.append((
            named_ancestors_reversed,
            set_value
        ))

    # Pop off all ancestors that are the same for all cell types.
    # e.g. 'cell', 'native cell', ...
    ancestor_list_lens = [len(x[0]) for x in ancestors_and_sets]
    min_ancestor_list_len = min(ancestor_list_lens)
    assert min_ancestor_list_len >= 1
    for level in range(min_ancestor_list_len - 1):
        unique_level_cell_types = set()
        for ancestors, cell_set in ancestors_and_sets:
            unique_level_cell_types.add(ancestors[0])

        if len(unique_level_cell_types) == 1:
            for ancestors, cell_set in ancestors_and_sets:
                ancestors.pop(0)
        else:
            break

    # Create the hierarchy as a dict.
    h = dict()
    for ancestors, cell_set in ancestors_and_sets:
        fill_in_dict_from_path(h, ancestors, cell_set)

    # Try removing all of the single-child "any" levels
    # now that the hierarchy has been created as a dict.
    h = remove_any_from_dict_levels(h)

    # Transform the dict into an object matching the JSON schema.
    tree["tree"] = [
        dict_to_tree("Cell Type Annotations", h)
    ]
    return tree
