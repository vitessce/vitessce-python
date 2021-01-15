import json
from uuid import uuid4

from .constants import (
    CoordinationType as ct,
    Component as cm,
    DataType as dt,
    FileType as ft
)

from .wrappers import (
    AnnDataWrapper,
    SnapWrapper,
)
from .widget import VitessceWidget
from .export import (
    export_to_s3,
    export_to_files,
)

def _get_next_scope(prev_scopes):
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    next_char_indices = [0]

    def next():
        r = []
        for char_index in next_char_indices:
            r = [chars[char_index]] + r
        increment = True
        for i in range(len(next_char_indices)):
            next_char_indices[i] += 1
            val = next_char_indices[i]
            if val >= len(chars):
                next_char_indices[i] = 0
            else:
                increment = False
                break
        
        if increment:
            next_char_indices.append(0)
        
        return "".join([ str(j) for j in r ])
    
    next_scope = next()
    while next_scope in prev_scopes:
        next_scope = next()
    
    return next_scope

class VitessceConfigDatasetFile:
    """
    A class to represent a file (described by a URL, data type, and file type) in a Vitessce view config dataset.
    """
    def __init__(self, url, data_type, file_type):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfigDataset.add_file()`` method.

        :param str url: A URL to this file. Can be a localhost URL or a remote URL.
        :param str data_type: A data type.
        :param str file_type: A file type.
        """
        self.file = {
            "url": url,
            "type": data_type,
            "fileType": file_type,
        }
    
    def to_dict(self):
        return self.file

class VitessceConfigDataset:
    """
    A class to represent a dataset (i.e. list of files containing common biological entities) in the Vitessce view config.
    """
    def __init__(self, uid, name):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfig.add_dataset()`` method.

        :param str uid: A unique identifier for this dataset.
        :param str name: A name for this dataset.
        """
        self.dataset = {
            "uid": uid,
            "name": name,
            "files": [],
        }
        self.objs = []
    
    def add_file(self, url, data_type, file_type):
        """
        Add a new file definition to this dataset instance.

        :param str url: The URL for the file, pointing to either a local or remote location.
        :param data_type: The type of data stored in the file. Must be compatible with the specified file type.
        :type data_type: str or vitessce.constants.DataType
        :param file_type: The file type. Must be compatible with the specified data type.
        :type file_type: str or vitessce.constants.FileType

        :returns: Self, to allow function chaining.
        :rtype: VitessceConfigDataset

        .. code-block:: python
            :emphasize-lines: 6-10

            from vitessce import VitessceConfig, DataType as dt, FileType as ft

            vc = VitessceConfig(name='My Config')
            my_dataset = (
                vc.add_dataset(name='My Dataset')
                .add_file(
                    url="http://example.com/cells.json",
                    data_type=dt.CELLS,
                    file_type=ft.CELLS_JSON,
                )
            )
        """

        assert type(data_type) == str or type(data_type) == dt
        assert type(file_type) == str or type(file_type) == ft

        # TODO: assert that the file type is compatible with the data type (and vice versa)

        if type(data_type) == str:
            data_type_str = data_type
        else:
            data_type_str = data_type.value
        
        if type(file_type) == str:
            file_type_str = file_type
        else:
            file_type_str = file_type.value

        self.dataset["files"].append(VitessceConfigDatasetFile(url=url, data_type=data_type_str, file_type=file_type_str))
        return self
    
    def add_object(self, obj):
        """
        Add a data object to this dataset instance.

        :param obj: A data object that can be served locally or uploaded to a remote storage provider.
        :type obj: anndata.AnnData or loompy.LoomConnection or zarr.Store

        :returns: Self, to allow function chaining.
        :rtype: VitessceConfigDataset
        """
        self.objs.append(obj)
        return self

    def to_dict(self, on_obj):
        obj_file_defs = []
        for obj_i, obj in enumerate(self.objs):
            if on_obj is not None:
                obj_file_defs += on_obj(obj, self.dataset["uid"], obj_i)
                
        return {
            **self.dataset,
            "files": [ f.to_dict() for f in self.dataset["files"] ] + obj_file_defs,
        }

class VitessceConfigViewHConcat:
    """
    A class to represent a horizontal concatenation of view instances.
    """
    def __init__(self, views):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``hconcat`` helper function.
        """
        self.views = views
    
    def __or__(self, other):
        return hconcat(self, other)
    
    def __truediv__(self, other):
        return vconcat(self, other)

def hconcat(*views):
    """
    A helper function to create a ``VitessceConfigViewHConcat`` instance.

    :param \*views: A variable number of views to concatenate horizontally.
    :type \*views: VitessceConfigView or VitessceConfigViewVConcat or VitessceConfigViewHConcat

    :returns: The concatenated view instance.
    :rtype: VitessceConfigViewHConcat

    .. code-block:: python
        :emphasize-lines: 8
        
        from vitessce import VitessceConfig, Component as cm, hconcat, vconcat

        vc = VitessceConfig()
        my_dataset = vc.add_dataset(name='My Dataset')
        v1 = vc.add_view(my_dataset, cm.SPATIAL)
        v2 = vc.add_view(my_dataset, cm.SPATIAL)
        v3 = vc.add_view(my_dataset, cm.SPATIAL)
        vc.layout(hconcat(v1, vconcat(v2, v3)))
    """
    return VitessceConfigViewHConcat(views)

class VitessceConfigViewVConcat:
    """
    A class to represent a vertical concatenation of view instances.
    """
    def __init__(self, views):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``vconcat`` helper function.
        """
        self.views = views
    
    def __or__(self, other):
        return hconcat(self, other)
    
    def __truediv__(self, other):
        return vconcat(self, other)

def vconcat(*views):
    """
    A helper function to create a ``VitessceConfigViewVConcat`` instance.

    :param \*views: A variable number of views to concatenate vertically.
    :type \*views: VitessceConfigView or VitessceConfigViewVConcat or VitessceConfigViewHConcat

    :returns: The concatenated view instance.
    :rtype: VitessceConfigViewVConcat

    .. code-block:: python
        :emphasize-lines: 8
        
        from vitessce import VitessceConfig, Component as cm, hconcat, vconcat

        vc = VitessceConfig()
        my_dataset = vc.add_dataset(name='My Dataset')
        v1 = vc.add_view(my_dataset, cm.SPATIAL)
        v2 = vc.add_view(my_dataset, cm.SPATIAL)
        v3 = vc.add_view(my_dataset, cm.SPATIAL)
        vc.layout(hconcat(v1, vconcat(v2, v3)))
    """
    return VitessceConfigViewVConcat(views)

class VitessceConfigView:
    """
    A class to represent a view (i.e. visualization component) in the Vitessce view config layout.
    """
    def __init__(self, component, coordination_scopes, x, y, w, h):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfig.add_view()`` method.

        :param str component: The name of the component used for this view.
        :param dict coordination_scopes: A mapping of coordination types to coordination scopes.
        :param int x: An x-coordinate for this view in the grid.
        :param int y: A y-coordinate for this view in the grid.
        :param int w: A width for this view in the grid.
        :param int h: A height for this view in the grid.
        """
        self.view = {
            "component": component,
            "coordinationScopes": coordination_scopes,
            "x": x,
            "y": y,
            "w": w,
            "h": h
        }
    
    def use_coordination(self, *c_scopes):
        """
        Attach a coordination scope to this view instance. All views using the same coordination scope for a particular coordination type will effectively be linked together.

        :param \*c_scopes: A variable number of coordination scope instances can be passed.
        :type \*c_scopes: VitessceConfigCoordinationScope

        :returns: Self, to allow chaining.
        :rtype: VitessceConfigView

        .. code-block:: python
            :emphasize-lines: 12-13

            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SPATIAL)
            zoom_scope, x_scope, y_scope = vc.add_coordination(
                ct.SPATIAL_ZOOM,
                ct.SPATIAL_TARGET_X,
                ct.SPATIAL_TARGET_Y,
            )
            v1.use_coordination(zoom_scope, x_scope, y_scope)
            v2.use_coordination(zoom_scope, x_scope, y_scope)
            zoom_scope.set_value(2)
            x_scope.set_value(0)
            y_scope.set_value(0)
        """
        for c_scope in c_scopes:
            assert type(c_scope) == VitessceConfigCoordinationScope
            self.view["coordinationScopes"][c_scope.c_type] = c_scope.c_scope
        return self
    
    def set_xywh(self, x, y, w, h):
        self.view["x"] = x
        self.view["y"] = y
        self.view["w"] = w
        self.view["h"] = h
        return self
    
    def set_props(self, **kwargs):
        if "props" in self.view.keys():
            self.view["props"] = {
                **self.view["props"],
                **kwargs
            }
        else:
            self.view["props"] = kwargs
        return self
    
    def to_dict(self):
        return self.view

    def __or__(self, other):
        return hconcat(self, other)
    
    def __truediv__(self, other):
        return vconcat(self, other)


class VitessceConfigCoordinationScope:
    """
    A class to represent a coordination scope in the Vitessce view config coordination space.
    """
    def __init__(self, c_type, c_scope):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfig.add_coordination()`` method.

        :param str c_type: The coordination type for this coordination scope.
        :param str c_scope: The coordination scope name.
        """
        self.c_type = c_type
        self.c_scope = c_scope
        self.c_value = None
    
    def set_value(self, c_value):
        """
        Set the value of the coordination scope.

        :param any c_value: The coordination value to be set. Can be any value that is valid for the coordination type. Must be serializable to JSON.

        :returns: Self, to allow chaining.
        :rtype: VitessceConfigCoordinationScope

        .. code-block:: python
            :emphasize-lines: 14-16

            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SPATIAL)
            zoom_scope, x_scope, y_scope = vc.add_coordination(
                ct.SPATIAL_ZOOM,
                ct.SPATIAL_TARGET_X,
                ct.SPATIAL_TARGET_Y,
            )
            v1.use_coordination(zoom_scope, x_scope, y_scope)
            v2.use_coordination(zoom_scope, x_scope, y_scope)
            zoom_scope.set_value(2)
            x_scope.set_value(0)
            y_scope.set_value(0)
        """
        self.c_value = c_value
        return self

class VitessceConfig:
    """
    A class to represent a Vitessce view config.
    """

    def __init__(self, name=None, description=None):
        """
        Construct a Vitessce view config object.

        :param str name: A name for the view config. Optional.
        :param str description: A description for the view config. Optional.

        .. code-block:: python
            :emphasize-lines: 3
            
            from vitessce import VitessceConfig

            vc = VitessceConfig(name='My Config')
        """
        self.config = {
            "version": "1.0.0",
            "name": name,
            "description": description,
            "datasets": [],
            "coordinationSpace": {},
            "layout": [],
            "initStrategy": "auto"
        }

        if name is None:
            self.config["name"] = ""
        else:
            self.config["name"] = name
        if description is None:
            self.config["description"] = ""
        else:
            self.config["description"] = description


    def add_dataset(self, name="", uid=None):
        """
        Add a dataset to the config.

        :param str name: A name for this dataset.
        :param str uid: A unique identifier for this dataset. Optional. If None, will be automatically generated.

        :returns: The instance for the new dataset.
        :rtype: VitessceConfigDataset
        
        .. code-block:: python
            :emphasize-lines: 5

            from vitessce import VitessceConfig, DataType as dt, FileType as ft

            vc = VitessceConfig(name='My Config')
            my_dataset = (
                vc.add_dataset(name='My Dataset')
                .add_file(
                    url="http://example.com/cells.json",
                    data_type=dt.CELLS,
                    file_type=ft.CELLS_JSON,
                )
            )
        """
        uid = uid if uid is not None else _get_next_scope([ d.dataset['uid'] for d in self.config["datasets"] ])
        assert type(uid) == str
        vcd = VitessceConfigDataset(uid, name)
        self.config["datasets"].append(vcd)
        [d_scope] = self.add_coordination(ct.DATASET)
        d_scope.set_value(uid)
        return vcd
    
    def add_view(self, dataset, component, x=0, y=0, w=1, h=1, mapping=None):
        """
        Add a view to the config.

        :param dataset: A dataset instance to be used for the data visualized in this view.
        :type dataset: VitessceConfigDataset
        :param component: A component name, either as a string or using the Component enum values.
        :type component: str or vitessce.constants.Component

        :param str mapping: An optional convenience parameter for setting the EMBEDDING_TYPE coordination scope value. This parameter is only applicable to the SCATTERPLOT component.
        :param int x: The horizontal position of the view. Must be an integer between 0 and 11. Optional.
        :param int y: The vertical position of the view. Must be an integer between 0 and 11. Optional.
        :param int w: The width of the view. Must be an integer between 1 and 12. Optional.
        :param int h: The height of the view. Must be an integer between 1 and 12. Optional.

        :returns: The instance for the new view.
        :rtype: VitessceConfigView

        .. code-block:: python
            :emphasize-lines: 5-6

            from vitessce import VitessceConfig, Component as cm

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SCATTERPLOT, mapping="X_umap")
        """
        assert type(dataset) == VitessceConfigDataset
        assert type(component) == str or type(component) == cm

        if type(component) == str:
            component_str = component
        else:
            component_str = component.value

        # Find the coordination scope name associated with the dataset
        dataset_matches = [
            scope_name
            for scope_name, dataset_scope in self.config["coordinationSpace"][ct.DATASET.value].items()
            if dataset_scope.c_value == dataset.dataset["uid"]
        ] if ct.DATASET.value in self.config["coordinationSpace"].keys() else []
        if len(dataset_matches) == 1:
            dataset_scope = dataset_matches[0]
        else:
            raise ValueError("No coordination scope matching the dataset parameter could be found in the coordination space.")

        # Set up the view's dataset coordination scope based on the dataset parameter.
        coordination_scopes = {
            ct.DATASET.value: dataset_scope,
        }
        vcv = VitessceConfigView(component_str, coordination_scopes, x, y, w, h)
        
        # Use the mapping parameter if component is scatterplot and the mapping is not None
        if mapping is not None:
            [et_scope] = self.add_coordination(ct.EMBEDDING_TYPE)
            et_scope.set_value(mapping)
            vcv.use_coordination(et_scope)
        self.config["layout"].append(vcv)
        return vcv
    

    def add_coordination(self, *c_types):
        """
        Add scope(s) for new coordination type(s) to the config.

        :param \*c_types: A variable number of coordination types.
        :type \*c_types: str or vitessce.constants.CoordinationType
        
        :returns: The instances for the new scope objects corresponding to each coordination type. These can be linked to views via the ``VitessceConfigView.use_coordination()`` method.
        :rtype: list[VitessceConfigCoordinationScope]

        .. code-block:: python
            :emphasize-lines: 7-11

            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SPATIAL)
            zoom_scope, x_scope, y_scope = vc.add_coordination(
                ct.SPATIAL_ZOOM,
                ct.SPATIAL_TARGET_X,
                ct.SPATIAL_TARGET_Y,
            )
            v1.use_coordination(zoom_scope, x_scope, y_scope)
            v2.use_coordination(zoom_scope, x_scope, y_scope)
            zoom_scope.set_value(2)
            x_scope.set_value(0)
            y_scope.set_value(0)
        """
        result = []
        for c_type in c_types:
            assert type(c_type) == ct or type(c_type) == str
            if type(c_type) == ct:
                c_type_str = c_type.value
            else:
                c_type_str = c_type
            prev_scopes = list(self.config["coordinationSpace"][c_type_str].keys()) if c_type_str in self.config["coordinationSpace"].keys() else []
            scope = VitessceConfigCoordinationScope(c_type_str, _get_next_scope(prev_scopes))
            if scope.c_type not in self.config["coordinationSpace"]:
                self.config["coordinationSpace"][scope.c_type] = {}
            self.config["coordinationSpace"][scope.c_type][scope.c_scope] = scope
            result.append(scope)
        return result
    
    def link_views(self, views, c_types, c_values = None):
        """
        A convenience function for setting up new coordination scopes across a set of views.
        
        :param views: views An array of view objects to link together.
        :type views: list of VitessceConfigView
        :param c_types: The coordination types on which to coordinate the views.
        :type c_types: list of str or list of vitessce.constants.CoordinationType
        :param list c_values: Initial values corresponding to each coordination type. Should have the same length as the c_types array. Optional.

        :returns: Self, to allow chaining.
        :rtype: VitessceConfig
        """
        c_scopes = self.add_coordination(*c_types)
        for view in views:
            for c_scope in c_scopes:
                view.use_coordination(c_scope)

        if c_values is not None and len(c_values) == len(c_types):
            for i, c_scope in enumerate(c_scopes):
                c_scope.set_value(c_values[i])

        return self
    
    def layout(self, view_concat):
        """
        Create a multi-view layout based on (potentially recursive) view concatenations.

        :param view_concat: Views arranged by concatenating vertically or horizontally. Alternatively, a single view can be passed.
        :type view_concat: VitessceConfigViewHConcat or VitessceConfigViewVConcat or VitessceConfigView

        :returns: Self, to allow chaining.
        :rtype: VitessceConfig
        
        .. code-block:: python
            :emphasize-lines: 8
            
            from vitessce import VitessceConfig, Component as cm, hconcat, vconcat

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SPATIAL)
            v3 = vc.add_view(my_dataset, cm.SPATIAL)
            vc.layout(hconcat(v1, vconcat(v2, v3)))
        
        .. code-block:: python
            :emphasize-lines: 8
            
            from vitessce import VitessceConfig, Component as cm

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            v2 = vc.add_view(my_dataset, cm.SPATIAL)
            v3 = vc.add_view(my_dataset, cm.SPATIAL)
            vc.layout(v1 | (v2 / v3)) # * magic * (alternative syntax)
        """

        def _layout(obj, x_min, x_max, y_min, y_max):
            w = x_max - x_min
            h = y_max - y_min
            if type(obj) == VitessceConfigView:
                obj.set_xywh(x_min, y_min, w, h)
            elif type(obj) == VitessceConfigViewHConcat:
                views = obj.views
                num_views = len(views)
                for i in range(num_views):
                    _layout(
                        views[i],
                        x_min+(w/num_views)*i,
                        x_min+(w/num_views)*(i+1),
                        y_min,
                        y_max
                    )
            elif type(obj) == VitessceConfigViewVConcat:
                views = obj.views
                num_views = len(views)
                for i in range(num_views):
                    _layout(
                        views[i],
                        x_min,
                        x_max,
                        y_min+(h/num_views)*i,
                        y_min+(h/num_views)*(i+1),
                    )

        # Recursively set the values (x,y,w,h) for each view.
        _layout(view_concat, 0, 12, 0, 12)

        # TODO: decide how to handle views that were omitted from the `view_concat` parameter
        # TODO: decide how to handle .add_view() being called after .layout() has been called

        return self
        
    def to_dict(self, on_obj=None):
        """
        Convert the view config instance to a dict object.

        :param on_obj: This callback is required only if datasets within the view config contain objects added via the ``VitessceConfigDataset.add_object`` method (rather than only files added via the ``VitessceConfigDataset.add_file`` method). This function must take the data object as a parameter, and return a list of valid file definition dicts (URL, data type, file type). This parameter is primarily intended to be used internally by the ``VitessceWidget`` class.
        :type on_obj: function or None

        :returns: The view config as a dict. Useful for serializing to JSON.
        :rtype: dict
        """
        return {
            **self.config,
            "datasets": [ d.to_dict(on_obj) for d in self.config["datasets"] ],
            "coordinationSpace": dict([
                (c_type, dict([
                    (c_scope_name, c_scope.c_value) for c_scope_name, c_scope in c_scopes.items() 
                ])) for c_type, c_scopes in self.config["coordinationSpace"].items()
            ]),
            # TODO: compute the x,y,w,h values if not explicitly defined
            "layout": [ c.to_dict() for c in self.config["layout"] ]
        }

    @staticmethod
    def from_dict(config):
        """
        Helper function to construct a Vitessce view config object from an existing config.

        :param dict config: An existing Vitessce view config as a dict to allow manipulation through the ``VitessceConfig`` API.
        
        :returns: The config instance.
        :rtype: VitessceConfig

        .. code-block:: python
            :emphasize-lines: 3
            
            from vitessce import VitessceConfig

            vc = VitessceConfig.from_dict(my_existing_config)
        """
        # TODO: Validate the incoming config.

        vc = VitessceConfig(name=config["name"], description=config["description"])

        # Add each dataset from the incoming config.
        for d in config["datasets"]:
            new_dataset = vc.add_dataset(uid=d["uid"], name=d["name"])
            for f in d["files"]:
                new_file = new_dataset.add_file(
                    url=f["url"],
                    data_type=f["type"],
                    file_type=f["fileType"]
                )
        
        for c_type in config['coordinationSpace'].keys():
            if c_type != ct.DATASET.value:
                c_obj = config['coordinationSpace'][c_type]
                vc.config['coordinationSpace'][c_type] = {}
                for c_scope_name, c_scope_value in c_obj.items():
                    scope = VitessceConfigCoordinationScope(c_type, c_scope_name)
                    scope.set_value(c_scope_value)
                    vc.config['coordinationSpace'][c_type][c_scope_name] = scope
        
        for c in config['layout']:
            c_coord_scopes = c['coordinationScopes'] if 'coordinationScopes' in c.keys() else {}
            new_view = VitessceConfigView(c['component'], c_coord_scopes, c['x'], c['y'], c['w'], c['h'])
            if 'props' in c.keys():
                new_view.set_props(**c['props'])
            vc.config['layout'].append(new_view)

        return vc


    @staticmethod
    def from_object(obj, name=None, description=None):
        """
        Helper function to automatically construct a Vitessce view config object from a single-cell dataset object.
        Particularly helpful when using the ``VitessceWidget`` Jupyter widget.

        :param obj: A single-cell dataset in the format of a commonly-used single-cell or imaging data analysis package.
        :type obj: anndata.AnnData or loompy.LoomConnection or zarr.Store
        
        :returns: The config instance.
        :rtype: VitessceConfig

        .. code-block:: python
            :emphasize-lines: 3
            
            from vitessce import VitessceConfig
            
            vc = VitessceConfig.from_object(my_scanpy_object)
        """
        vc = VitessceConfig(name=name, description=description)
        dataset = vc.add_dataset()

        if type(obj) == AnnDataWrapper:
            dataset.add_object(obj)

            # TODO: use the available embeddings to determine how many / which scatterplots to add.
            scatterplot = vc.add_view(dataset, cm.SCATTERPLOT, mapping="X_umap")
            cell_sets = vc.add_view(dataset, cm.CELL_SETS)
            genes = vc.add_view(dataset, cm.GENES)
            heatmap = vc.add_view(dataset, cm.HEATMAP)

            vc.layout((scatterplot | (cell_sets / genes)) / heatmap)
        
        elif type(obj) == SnapWrapper:
            dataset.add_object(obj)

            genomic_profiles = vc.add_view(dataset, cm.GENOMIC_PROFILES)
            scatter = vc.add_view(dataset, cm.SCATTERPLOT, mapping = "UMAP")
            cell_sets = vc.add_view(dataset, cm.CELL_SETS)

            vc.layout(genomic_profiles / (scatter | cell_sets))
        
        else:
            print("Encountered an unknown object type. Unable to automatically generate a Vitessce view config.")

        return vc
    
    def widget(self, **kwargs):
        """
        Convenience function for instantiating a VitessceWidget object based on this config.
        
        :param str theme: The theme name, either "light" or "dark". By default, "auto", which selects light or dark based on operating system preferences.
        :param int height: The height of the widget, in pixels. By default, 600.
        :param int port: The port to use when serving data objects on localhost. By default, 8000.
        :param bool proxy: Is this widget being served through a proxy, for example with a cloud notebook (e.g. Binder)?
        
        :returns: The Jupyter widget.
        :rtype: VitessceWidget

        .. code-block:: python
            :emphasize-lines: 6-7

            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            vc.layout(v1)
            vw = vc.widget()
            vw
        """
        return VitessceWidget(self, **kwargs)
        
    def export(self, to, *args, **kwargs):
        """
        Export this config's data objects to the local file system or a cloud storage system and get the resulting view config.
        
        :param str to: The export destination. Valid values include "S3" and "files".
        :param \*\*kwargs: Keyword arguments to pass to the export function.
        :returns: The config as a dict, with URLs for the bucket filled in.
        :rtype: dict
        
        .. code-block:: python
            :emphasize-lines: 8
        
            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(my_dataset, cm.SPATIAL)
            vc.layout(v1)
            
            config_dict = vc.export(to="S3")
        """
        if to == "S3":
            return export_to_s3(self, *args, **kwargs)
        elif to == "files":
            return export_to_files(self, *args, **kwargs)
        else:
            raise ValueError("Unknown export destination.")


