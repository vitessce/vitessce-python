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
	try:
		cell_ids = adata.obs.index.tolist()
		cluster_ids = adata.obs['ClusterID'].unique().tolist()
		cell_cluster_ids = adata.obs['ClusterID'].values.tolist()
		
		cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))
		
		cell_sets_json = {
			"datatype": "cell",
			"version": "0.1.2",
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
	except Exception as e:
		print(str(e))
			
	#print(cell_sets_json)
	
	async def response_func(req):
		try:
			return JSONResponse(cell_sets_json, headers=HEADERS)
		except exception as e:
			print(str(e))
			return JSONResponse({}, headers=HEADERS)
		
	return response_func

def create_config_and_routes(config):

	routes = []

	def create_obj_route(obj):
		dtype = type(obj)
		print(f"Creating server route for object with type {dtype}")
		
		# Check if dtype is AnnData
		try:
			import anndata
			if dtype == anndata.AnnData:
				routes.append(Route('/my-dataset/cells', create_anndata_response_cells_json(data)))
				routes.append(Route('/my-dataset/cell-sets', create_anndata_response_cell_sets_json(data)))
				return [
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
		except ImportError:
			pass
		
		# Check if dtype is Loom
		try:
			import loompy
			if dtype == loompy.LoomConnection:
				# TODO: append routes
				return [
					# TODO: add file definitions
				]
		except ImportError:
			pass
	
	config_dict = config.to_dict(on_obj=create_obj_route)
	return config_dict, routes