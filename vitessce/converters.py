class Cells:

  def __init__(self, cell_ids=[]):
    self.cell_ids = cell_ids
    self.json = dict(zip(cell_ids, [{} for _ in cell_ids]))

  def add_mapping(self, name, coords):
    for idx, id in enumerate(self.cell_ids):
      if 'mappings' not in self.json[id]:
        self.json[id]['mappings'] = { name: coords[idx] }
      else:
        self.json[id]['mappings'][name] = coords[idx]
  
  def add_genes(self, genes):
    for idx, id in enumerate(self.cell_ids):
      self.json[id]['genes'] = genes[idx]
  
  def add_xy(self, xy):
    for idx, id in enumerate(self.cell_ids):
      self.json[id]['xy'] = xy[idx]
  
  def add_factors(self, factors):
    for idx, id in enumerate(self.cell_ids):
      self.json[id]['factors'] = factors[idx]

  def add_polygon_outline(self, polygon_outline):
    for idx, id in enumerate(self.cell_ids):
      self.json[id]['poly'] = polygon_outline[idx]

class CellSets:

  def __init__(self, first_node_name):
    self.json = {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": [{
            "name": first_node_name,
            "children": []
        }]
    }

  def add_set(self, name, parent_path, cell_set = None):
    parent_node = self.tree_find_node_by_path(parent_path)
    new_node = { "name": name }
    if cell_set:
      new_node['set'] = cell_set
    if 'children' not in parent_node:
      parent_node['children'] = [new_node]
    else:
      parent_node['children'].append(new_node)
  
  def find_node_by_path(self, node, path, curr_index):
    curr_node_name = path[curr_index]
    if node['name'] == curr_node_name:
      if curr_index == len(path) - 1:
        return node
      if 'children' in node:
        found_nodes = [
          self.find_node_by_path(child, path, curr_index + 1) for child in node['children']
        ]
        found_nodes_not_none = [n for n in found_nodes if n]
        if len(found_nodes_not_none) == 1:
          return found_nodes[0]
    return None

  def tree_find_node_by_path(self, path):
    found_nodes = [self.find_node_by_path(node, path, 0) for node in self.json['tree']]
    found_nodes_not_none = [n for n in found_nodes if n]
    if len(found_nodes_not_none) == 1:
      return found_nodes_not_none[0]
    return None