from starlette.responses import JSONResponse, UJSONResponse
from starlette.routing import Route

HEADERS = { 'Access-Control-Allow-Origin': '*' }

def create_dummy_response_cells_json():

	async def response_func(req):
		return JSONResponse({"hello": "cells"}, headers=HEADERS)
	
	print("done creating response")
		
	return response_func

def create_dummy_response_cell_sets_json():

	async def response_func(req):
		return JSONResponse({"hello": "cell sets"}, headers=HEADERS)
	
	print("done creating response")
		
	return response_func


def create_anndata_json(adata):
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

	cell_sets_json = {
		"datatype": "cell",
		"version": "0.1.2",
		"tree": [{
			"name": "Clusters",
			"children": []
		}]
	}

	cell_ids = adata.obs.index.tolist()
	cluster_ids = adata.obs['CellType'].unique().tolist()
	cell_cluster_ids = adata.obs['CellType'].values.tolist()
	
	cell_cluster_tuples = list(zip(cell_ids, cell_cluster_ids))
		
	for cluster_id in cluster_ids:
		cell_sets_json["tree"][0]["children"].append({
			"name": str(cluster_id),
			"set": [
				str(cell_id)
				for cell_id, cell_cluster_id in cell_cluster_tuples
				if cell_cluster_id == cluster_id
			]
		})

	return cells_json, cell_sets_json

def create_response_json(data_json):
	async def response_func(req):
		return UJSONResponse(data_json, headers=HEADERS)
	return response_func


def create_obj_routes(obj, obj_i, dataset_uid):
	"""
	For a particular data object, simultaneously set up:
	
	* its server routes and their responses, and
	* the corresponding view config dataset file definitions

	:param obj: An object representing a single-cell data analysis result or microscopy image.
	:type obj: anndata.AnnData or loompy.Loom or zarr.Store
	
	:returns: A list of view config file definitions and a list of server routes.
	:rtype: tuple[list[dict], list[starlette.routing.Route]]
	"""
	obj_file_defs = []
	obj_routes = []

	dtype = type(obj)
	print(f"Creating server route for object with type {dtype}")
	
	# Check if dtype is AnnData
	try:
		import anndata
		if dtype == anndata.AnnData:
			cells_json, cell_sets_json = create_anndata_json(obj)

			obj_routes = [
				Route(f'/{dataset_uid}/{obj_i}/cells', create_response_json(cells_json)),
				Route(f'/{dataset_uid}/{obj_i}/cell-sets', create_response_json(cell_sets_json)),
			]
			obj_file_defs = [
				{
					"type": "cells",
					"fileType": "cells.json",
					"url": f"http://localhost:8000/{dataset_uid}/{obj_i}/cells"
				},
				{
					"type": "cell-sets",
					"fileType": "cells-sets.json",
					"url": f"http://localhost:8000/{dataset_uid}/{obj_i}/cell-sets"
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

	
	return obj_file_defs, obj_routes

def create_exception_handlers():

	async def not_found(request, exc):
		return JSONResponse({"error": "not found"}, status_code=exc.status_code, headers=HEADERS)

	async def server_error(request, exc):
		return JSONResponse({"error": "server error"}, status_code=exc.status_code, headers=HEADERS)


	exception_handlers = {
		404: not_found,
		500: server_error
	}
	return exception_handlers

def create_dummy_routes():
	return [
		Route('/my-dataset/cells', create_dummy_response_cells_json()),
		Route('/my-dataset/cell-sets', create_dummy_response_cell_sets_json())
	]
			
	
