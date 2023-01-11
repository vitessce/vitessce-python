import networkx
import obonet

from constants import CL_ROOT_ID


def load_cl_obo_graph(cl_obo_file):
    graph = obonet.read_obo(cl_obo_file)

    # Make sure there are no cycles.
    assert networkx.is_directed_acyclic_graph(graph)

    id_to_name = {
        id_: data.get('name')
        for id_, data in graph.nodes(data=True)
    }
    name_to_id = {
        data['name']: id_
        for id_, data in graph.nodes(data=True) if 'name' in data
    }

    return graph, id_to_name, name_to_id


# Construct the tree, according to the following schema:
# https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
def init_cell_sets_tree():
    return {
        "datatype": "obs",
        "version": "0.1.3",
        "tree": []
    }


# Merge multiple cell sets trees.
def merge_cell_sets_trees(*args):
    '''
    >>> tree_a = init_cell_sets_tree()
    >>> tree_b = init_cell_sets_tree()
    >>> tree_a["tree"].append({
    ...   "name": "a",
    ...   "set": [1, 2, 3]
    ... })
    >>> tree_b["tree"].append({
    ...   "name": "b",
    ...   "set": [4, 5, 6]
    ... })
    >>> tree_merged = merge_cell_sets_trees(tree_a, tree_b)
    >>> len(tree_merged["tree"])
    2
    >>> [x["name"] for x in tree_merged["tree"]]
    ['a', 'b']
    '''
    assert len(args) > 1

    result = args[0]
    for arg in args[1:]:
        for level_zero_node in arg["tree"]:
            result["tree"].append(level_zero_node)
    return result


# Recursively convert a nested dict to a level zero node
# of the cell-set hierarchy schema.
def dict_to_tree(name, value):
    '''
    >>> h_dict = {
    ...     "hematopoietic cell": {
    ...         "leukocyte": [4, 5, 6],
    ...         "hematopoietic precursor cell": [7, 8, 9]
    ...     },
    ...     "epithelial cell": [1, 2, 3]
    ... }
    >>> h_tree = dict_to_tree("test", h_dict)
    >>> h_tree["name"]
    'test'
    >>> [x["name"] for x in h_tree["children"]]
    ['hematopoietic cell', 'epithelial cell']
    '''
    if isinstance(value, dict):
        return {
            "name": name,
            "children": [
                dict_to_tree(child_name, child_value)
                for child_name, child_value in value.items()
            ]
        }
    else:
        return {
            "name": name,
            "set": value,
        }


# Given a list of multiple paths up the DAG,
# sort the list according to a heuristic.
def sort_paths_up_cell_ontology(paths_up):
    '''
    >>> ex_paths_up = [
    ...     ['b', 'motile cell', 'native cell', 'cell'],
    ...     ['a', 'b', 'c', 'somatic cell', 'native cell', 'cell'],
    ...     ['a', 'somatic cell', 'native cell', 'cell']
    ... ]
    >>> sorted_paths_up = sort_paths_up_cell_ontology(ex_paths_up)
    >>> sorted_paths_up[0]
    ['a', 'somatic cell', 'native cell', 'cell']
    >>> sorted_paths_up[1]
    ['a', 'b', 'c', 'somatic cell', 'native cell', 'cell']
    '''
    PREFERENCES = [
        ['animal cell', 'eukaryotic cell', 'native cell', 'cell'],
        ['somatic cell', 'native cell', 'cell'],
        ['nucleate cell', 'native cell', 'cell'],
        ['precursor cell', 'native cell', 'cell'],
    ]
    # Prefer all of the above before "functional" categories like
    # [..., 'motile cell', 'native cell', 'cell']
    WORST_PREFERENCE_INDEX = len(PREFERENCES)

    def get_first_preference_index_and_path_length(path_up):
        path_preference_match_index = WORST_PREFERENCE_INDEX
        for preference_index, preference in enumerate(PREFERENCES):
            if path_up[-len(preference):] == preference:
                path_preference_match_index = preference_index
                break
        # Return a tuple of the first matching preference "path ending"
        # and the path length (to use shorter paths if multiple paths
        # match the same top path ending).
        return (path_preference_match_index, len(path_up))
    return sorted(paths_up, key=get_first_preference_index_and_path_length)


# Using the cell ontology DAG,
# get all possible paths up to the root from
# a given node_id value.
def get_paths_up_cell_ontology(graph, node_id):
    if node_id == CL_ROOT_ID:
        return [
            [node_id]
        ]

    # Get ancestors of the cell type
    # (counterintuitive that the function is called descendants).
    ancestor_term_set = networkx.descendants(graph, node_id)

    # Make sure the cell type has an ancestor
    # with the 'cell' root ID.
    assert CL_ROOT_ID in ancestor_term_set

    # Get the parents of the current node.
    node_parents = graph.out_edges(node_id, keys=True)

    up_dag_paths = []
    for node_parent in node_parents:
        _, curr_parent_id, relationship = node_parent
        if relationship == "is_a":
            parent_paths = get_paths_up_cell_ontology(graph, curr_parent_id)
            for parent_path in parent_paths:
                up_dag_paths.append([node_id] + parent_path)
    return up_dag_paths


# Using a path through the DAG for a particular cell set,
# recursively fill in a nested dict.
def fill_in_dict_from_path(d, keys, child):
    '''
    >>> h = dict()
    >>> path_set_tuple = (
    ...     ["type a", "type b", "type c"],
    ...     ["cell 1", "cell 2", "cell 3"]
    ... )
    >>> fill_in_dict_from_path(h, path_set_tuple[0], path_set_tuple[1])
    {'any': ['cell 1', 'cell 2', 'cell 3']}
    >>> h
    {'type a': {'type b': {'type c': {'any': ['cell 1', 'cell 2', 'cell 3']}}}}
    '''
    """
    d : dict The resulting dictionary so far.
    keys : A list of keys representing the path down the cell ontology DAG.
    child : A set value.
    """
    key = keys[0]

    if key in d and isinstance(d[key], dict):
        result = d[key]
    else:
        result = d[key] = dict()

    if len(keys) == 1:
        result["any"] = child
        return result
    else:
        new_keys = keys.copy()
        new_keys.pop(0)
        return fill_in_dict_from_path(result, new_keys, child)


# Try removing the extra hierarchy level for the "any" set,
# if it is an "only-child"
def remove_any_from_dict_levels(v):
    '''
    >>> h = {'type b': {'type c': {'any': ['cell 1', 'cell 2', 'cell 3']}}}
    >>> new_h = remove_any_from_dict_levels(h)
    >>> new_h
    {'type b': {'type c': ['cell 1', 'cell 2', 'cell 3']}}
    '''
    if isinstance(v, dict):
        keys = list(v.keys())
        if len(keys) == 1 and keys[0] == "any":
            # Return the value associated with the "any" property,
            # since "any" has no siblings.
            return v["any"]
        else:
            # This is a dict with multiple values, so recursively
            # try this function on all of its values.
            return dict(
                zip(
                    keys,
                    list(map(remove_any_from_dict_levels, v.values()))
                )
            )
    else:
        # This is not a dict, so just return as-is.
        return v
