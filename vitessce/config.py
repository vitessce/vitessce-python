from starlette.responses import JSONResponse
from starlette.routing import Route

HEADERS = { 'Access-Control-Allow-Origin': '*' }

def create_anndata_response_cells_json(adata):
	available_embeddings = list(adata.obsm.keys())
	
	cell_ids = adata.obs.index.tolist()
	cell_mappings = []
	for e in available_embeddings:
		mapping = adata.obsm[e][:,0:2].tolist()
		cell_mappings.append(list(zip(
			[e for i in range(len(mapping))],
			mapping
		)))
	cell_mappings_zip = list(zip(*cell_mappings))
	cells_json = dict(zip(
		cell_ids,
		[
			{ 'mappings': dict(cell_mapping), 'genes': {} }
			for cell_mapping in cell_mappings_zip
		]
	))
	
	async def response_func(req):
		return JSONResponse(cells_json, headers=HEADERS)
	
	return response_func

def create_anndata_response_cell_sets_json(adata):
	cell_ids = adata.obs.index.tolist()
	cluster_ids = adata.obs['ClusterID'].unique().tolist()
	cell_cluster_ids = adata.obs['ClusterID'].values.tolist()
	
	cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))
	
	cell_sets_json = {
		"datatype": "cell",
		"version": "0.1.3",
		"tree": [{
			"name": "Clusters",
			"children": []
		}]
	}
		
	for cluster_id in cluster_ids:
		cell_sets_json["tree"][0]["children"].append({
			"name": str(cluster_id),
			"set": [
				str(cell_id)
				for cell_id, cell_cluster_id in cell_cluster_tuples
				if cell_cluster_id == cluster_id
			]
		})
	
	async def response_func(req):
		return JSONResponse(cell_sets_json, headers=HEADERS)
	
	return response_func

def create_config_and_routes(data):
	dtype = type(data)
	print(dtype)
	
	config = {
		"version": "1.0.0",
		"name": "",
		"description": "",
		"datasets": [],
		"coordinationSpace": {
			"embeddingType": {}
		},
		"layout": [],
		"initStrategy": "auto"
		
	}
	routes = []
	
	if dtype == list:
		for item in data:
			item_config, item_routes = create_config_and_routes(item)
			# TODO: merge multiple configs / routes into single coordinated config and list of routes
	
	else:
		# Check if dtype is AnnData
		try:
			import anndata
			if dtype == anndata.AnnData:
				available_embeddings = list(data.obsm.keys())
				
				config["datasets"].append(
					{
						"uid": "my-dataset",
						"files": [
							{
								"type": "cells",
								"fileType": "cells.json",
								"url": "http://localhost:8000/my-dataset/cells"
							},
							{
								"type": "cell-sets",
								"fileType": "cells-sets.json",
								"url": "http://localhost:8000/my-dataset/cell-sets"
							}
						]
					}
				)
				for i, e in enumerate(available_embeddings):
					config["coordinationSpace"]["embeddingType"][e] = e
					config["layout"].append({
						"component": 'scatterplot',
						"coordinationScopes": {
							"embeddingType": e,
						},
						"x": i*4, "y": 0, "w": 4, "h": 2
					})
				routes.append(Route('/my-dataset/cells', create_anndata_response_cells_json(data)))
				config["layout"].append({
					"component": "cellSets",
					"coordinationScopes": {
						
					},
					"x": 0, "y": 2, "w": 4, "h": 2
				})
				routes.append(Route('/my-dataset/cell-sets', create_anndata_response_cell_sets_json(data)))
		except ImportError:
			pass
		
		# Check if dtype is Loom
		try:
			import loompy
			if dtype == loompy.LoomConnection:
				pass # TODO: support Loom files
		except ImportError:
			pass
	
	
	return config, routes