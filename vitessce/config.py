from starlette.responses import JSONResponse
from starlette.routing import Route
import anndata
import loompy

HEADERS = { 'Access-Control-Allow-Origin': '*' }

async def homepage(req):
	return JSONResponse({'hello': 'world'}, headers=HEADERS)

def create_config_and_routes(data):
	dtype = type(data)
	print(dtype)
	
	config = {
		"version": "1.0.0",
		"name": "COVID-19 patient PBMCs",
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
	
	elif dtype == anndata.AnnData:
		available_embeddings = list(data.obsm.keys())
		
		config["datasets"].append(
			{
				"uid": "my-dataset",
				"files": [
					{
						"type": "cells",
						"fileType": "cells.json",
						"url": "http://localhost:8000/my-dataset/cells"
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
				"x": i*4, "y": 0, "w": 3, "h": 2
			})
		routes.append(Route('/my-dataset/cells', homepage))
	
	elif dtype == loompy.LoomConnection:
		pass # TODO: support Loom files
	
	
	return config, routes