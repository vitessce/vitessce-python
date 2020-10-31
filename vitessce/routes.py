
from .constants import DataType as dt

def create_obj_routes(obj, port, dataset_uid, obj_i):
	"""
	For a particular data object, simultaneously set up:
	
	* its server routes and their responses, and
	* the corresponding view config dataset file definitions

	:param obj: An object representing a single-cell data analysis result or microscopy image.
	:type obj: anndata.AnnData or loompy.LoomConnection or zarr.hierarchy.Group
	
	:returns: A list of view config file definitions and a list of server routes.
	:rtype: tuple[list[dict], list[starlette.routing.Route]]
	"""
	obj_file_defs = []
	obj_routes = []

	for data_type in dt:
		try:
			dt_obj_file_defs, dt_obj_routes = obj._get_data(data_type, port, dataset_uid, obj_i)
			obj_file_defs += dt_obj_file_defs
			obj_routes += dt_obj_routes
		except NotImplementedError:
			pass

	return obj_file_defs, obj_routes
	