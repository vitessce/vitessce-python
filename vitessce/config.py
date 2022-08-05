import sys
import inspect
import copy as copy_module
import black
from collections import OrderedDict

from .constants import (
    CoordinationType as ct,
    Component as cm,
    DataType as dt,
    FileType as ft
)

from .repr import make_repr, make_params_repr


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

        return "".join([str(j) for j in r])

    next_scope = next()
    while next_scope in prev_scopes:
        next_scope = next()

    return next_scope


class VitessceConfigDatasetFile:
    """
    A class to represent a file (described by a URL, data type, and file type) in a Vitessce view config dataset.
    """

    def __init__(self, data_type, file_type, url=None, options=None):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfigDataset.add_file()`` method.

        :param str data_type: A data type.
        :param str file_type: A file type.
        :param url: A URL to this file. Can be a localhost URL or a remote URL.
        :type url: str or None
        :param options: Extra options to pass to the file loader class.
        :type options: dict or list or None
        """
        self.file = {
            "type": data_type,
            "fileType": file_type
        }
        if url:
            self.file["url"] = url
        if options:
            self.file["options"] = options

    def __repr__(self):
        repr_dict = {
            "data_type": self.file["type"],
            "file_type": self.file["fileType"],
        }
        if "url" in self.file:
            repr_dict["url"] = self.file["url"]
        if "options" in self.file:
            repr_dict["options"] = self.file["options"]

        return make_repr(repr_dict, class_def=self.__class__)

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

    def _to_py_params(self):
        return {
            "uid": self.dataset["uid"],
            "name": self.dataset["name"]
        }

    def get_name(self):
        """
        Get the name for this dataset.

        :returns: The name.
        :rtype: str
        """
        return self.dataset["name"]

    def get_uid(self):
        """
        Get the uid value for this dataset.

        :returns: The uid.
        :rtype: str
        """
        return self.dataset["uid"]

    def add_file(self, data_type, file_type, url=None, options=None):
        """
        Add a new file definition to this dataset instance.

        :param data_type: The type of data stored in the file. Must be compatible with the specified file type.
        :type data_type: str or vitessce.constants.DataType
        :param file_type: The file type. Must be compatible with the specified data type.
        :type file_type: str or vitessce.constants.FileType
        :param url: The URL for the file, pointing to either a local or remote location.
        :type url: str or None
        :type options: Extra options to pass to the file loader class. Optional.
        :type options: dict or list or None

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

        assert isinstance(data_type, str) or isinstance(data_type, dt)
        assert isinstance(file_type, str) or isinstance(file_type, ft)

        # TODO: assert that the file type is compatible with the data type (and vice versa)

        if isinstance(data_type, str):
            data_type_str = data_type
        else:
            data_type_str = data_type.value

        if isinstance(file_type, str):
            file_type_str = file_type
        else:
            file_type_str = file_type.value

        self._add_file(VitessceConfigDatasetFile(
            url=url, data_type=data_type_str, file_type=file_type_str, options=options))
        return self

    def _add_file(self, obj):
        self.dataset["files"].append(obj)
        return self

    def add_object(self, obj):
        """
        Add a data object to this dataset instance.

        :param obj: A data object that can be served locally or uploaded to a remote storage provider.
        :type obj: anndata.AnnData or loompy.LoomConnection or zarr.Store

        :returns: Self, to allow function chaining.
        :rtype: VitessceConfigDataset
        """
        obj.convert_and_save(self.dataset["uid"], len(self.objs))
        self.objs.append(obj)
        return self

    def _get_files(self):
        """
        Get a list of files and data objects associated with this dataset.

        :returns: The list of files and datasets.
        :rtype: list of VitessceConfigDatasetFile
        """
        return self.dataset["files"]

    def _get_objects(self):
        """
        Get a list of data objects associated with this dataset.

        :returns: The list of data objects.
        :rtype: list of AbstractWrapper instances
        """
        return self.objs

    def to_dict(self, base_url=None):
        obj_file_defs = []
        for obj in self.objs:
            obj_file_defs += obj.get_file_defs(base_url)

        return {
            **self.dataset,
            "files": [f.to_dict() for f in self.dataset["files"]] + obj_file_defs,
        }

    def get_routes(self):
        routes = []
        for obj in self.objs:
            routes += obj.get_routes()

        return routes


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
        v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
        v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
        v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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
        v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
        v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
        v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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

    def _to_py_params(self):
        params_dict = {
            "component": self.view["component"],
            "x": self.view["x"],
            "y": self.view["y"],
            "w": self.view["w"],
            "h": self.view["h"]
        }
        # Only include coordination_scopes if there are coordination scopes other than
        # the coorindation scope for the 'dataset' coordination type.
        non_dataset_coordination_scopes = {
            c_type: c_scope
            for c_type, c_scope in self.view["coordinationScopes"].items()
            if c_type != ct.DATASET.value
        }
        if len(non_dataset_coordination_scopes) > 0:
            params_dict["coordination_scopes"] = non_dataset_coordination_scopes
        return params_dict

    def get_coordination_scope(self, c_type):
        """
        Get the coordination scope name for a particular coordination type.

        :param str c_type: The coordination type of interest.

        :returns: The coordination scope name.
        :rtype: str or None
        """
        return self.view["coordinationScopes"].get(c_type)

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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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
        """
        Set the dimensions for this view.

        :param int x: The horizontal position.
        :param int y: The vertical position.
        :param int w: The width.
        :param int h: The height.

        :returns: Self, to allow chaining.
        :rtype: VitessceConfigView
        """
        self.view["x"] = x
        self.view["y"] = y
        self.view["w"] = w
        self.view["h"] = h
        return self

    def set_props(self, **kwargs):
        """
        Set the props for this view.

        :param \*\*kwargs: A variable number of named props.

        :returns: Self, to allow chaining.
        :rtype: VitessceConfigView
        """
        if "props" in self.view:
            self.view["props"].update(kwargs)
        else:
            self.view["props"] = kwargs
        return self

    def get_props(self):
        """
        Get the props for this view.

        :returns: The props.
        :rtype: dict or None
        """
        return self.view.get("props")

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

    def __init__(self, c_type, c_scope, c_value=None):
        """
        Not meant to be instantiated directly, but instead created and returned by the ``VitessceConfig.add_coordination()`` method.

        :param str c_type: The coordination type for this coordination scope.
        :param str c_scope: The coordination scope name.
        :param c_value: The value for the coordination scope. Optional.
        """
        self.c_type = c_type
        self.c_scope = c_scope
        self.c_value = c_value

    def _to_py_params(self):
        return {
            "c_type": self.c_type,
            "c_scope": self.c_scope,
            "c_value": self.c_value,
        }

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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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

    def __init__(self, name=None, description=None, schema_version="1.0.7"):
        """
        Construct a Vitessce view config object.

        :param str name: A name for the view config. Optional.
        :param str description: A description for the view config. Optional.
        :param str schema_version: The view config schema version.

        .. code-block:: python
            :emphasize-lines: 3

            from vitessce import VitessceConfig

            vc = VitessceConfig(name='My Config')
        """
        self.config = {
            "version": schema_version,
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

    def _to_py_params(self):
        return {
            "name": self.config["name"],
            "description": self.config["description"],
            "schema_version": self.config["version"],
        }

    def add_dataset(self, name="", uid=None, files=None, objs=None):
        """
        Add a dataset to the config.

        :param str name: A name for this dataset.
        :param str uid: A unique identifier for this dataset. Optional. If None, one will be automatically generated.
        :param files: A list of VitessceConfigDatasetFile instances. optional.
        :type files: list or None
        :param objs: A list of AbstractWrapper instances. Optional.
        :type objs: list or None

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
        uid = uid if uid is not None else _get_next_scope(
            [d.dataset['uid'] for d in self.config["datasets"]])
        assert type(uid) == str
        vcd = VitessceConfigDataset(uid, name)
        self.config["datasets"].append(vcd)
        [d_scope] = self.add_coordination(ct.DATASET)
        d_scope.set_value(uid)

        if isinstance(files, list):
            for obj in files:
                vcd._add_file(obj)
        if isinstance(objs, list):
            for obj in objs:
                vcd.add_object(obj)

        return vcd

    def get_dataset_by_uid(self, uid):
        """
        Get a dataset associated with this configuration based on its uid.

        :param str uid: The unique identifier for the dataset of interest.

        :returns: The dataset object.
        :rtype: VitessceConfigDataset or None
        """
        for dataset in self.config["datasets"]:
            if dataset.get_uid() == uid:
                return dataset
        return None

    def get_dataset_by_coordination_scope_name(self, query_scope_name):
        """
        Get a dataset associated with this configuration based on a coordination scope.

        :param str query_scope_name: The unique identifier for the dataset coordination scope of interest.

        :returns: The dataset object.
        :rtype: VitessceConfigDataset or None
        """
        if ct.DATASET.value in self.config["coordinationSpace"]:
            for scope_name, dataset_scope in self.config["coordinationSpace"][ct.DATASET.value].items():
                if scope_name == query_scope_name:
                    return self.get_dataset_by_uid(dataset_scope.c_value)
        return None

    def get_datasets(self):
        """
        Get the datasets associated with this configuration.

        :returns: The list of dataset objects.
        :rtype: list of VitessceConfigDataset
        """
        return self.config["datasets"]

    def add_view(self, component, dataset=None, dataset_uid=None, x=0, y=0, w=1, h=1, mapping=None, coordination_scopes=None, props=None):
        """
        Add a view to the config.

        :param component: A component name, either as a string or using the Component enum values.
        :type component: str or vitessce.constants.Component
        :param dataset: A dataset instance to be used for the data visualized in this view. Must provide dataset or dataset_uid, but not both.
        :type dataset: VitessceConfigDataset or None
        :param dataset_uid: A unique ID for a dataset to be used for the data visualized in this view. Must provide dataset or dataset_uid, but not both.
        :type dataset_uid: str or None

        :param str mapping: An optional convenience parameter for setting the EMBEDDING_TYPE coordination scope value. This parameter is only applicable to the SCATTERPLOT component.
        :param int x: The horizontal position of the view. Must be an integer between 0 and 11. Optional. This will be ignored if you call the `layout` method of this class using `VitessceConfigViewHConcat` and `VitessceConfigViewVConcat` objects.
        :param int y: The vertical position of the view. Must be an integer between 0 and 11. Optional. This will be ignored if you call the `layout` method of this class using `VitessceConfigViewHConcat` and `VitessceConfigViewVConcat` objects.
        :param int w: The width of the view. Must be an integer between 1 and 12. Optional. This will be ignored if you call the `layout` method of this class using `VitessceConfigViewHConcat` and `VitessceConfigViewVConcat` objects.
        :param int h: The height of the view. Must be an integer between 1 and 12. Optional. This will be ignored if you call the `layout` method of this class using `VitessceConfigViewHConcat` and `VitessceConfigViewVConcat` objects.
        :param coordination_scopes: A mapping from coordination types to coordination scope names for this view.
        :type coordination_scopes: dict or None
        :param props: Props to set for the view using the VitessceConfigView.set_props method.
        :type props: dict or None

        :returns: The instance for the new view.
        :rtype: VitessceConfigView

        .. code-block:: python
            :emphasize-lines: 5-6

            from vitessce import VitessceConfig, Component as cm

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SCATTERPLOT, dataset=my_dataset, mapping="X_umap")
        """
        # User should only provide dataset or dataset_uid, but not both.
        assert isinstance(dataset, VitessceConfigDataset) or isinstance(
            dataset_uid, str)
        assert dataset is None or dataset_uid is None
        assert type(component) in [str, cm]

        if dataset is None:
            dataset = self.get_dataset_by_uid(dataset_uid)
            if dataset is None:
                raise ValueError(
                    "A dataset with the provided dataset_uid could not be found.")

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
            raise ValueError(
                "No coordination scope matching the dataset parameter could be found in the coordination space.")

        # Set up the view's dataset coordination scope based on the dataset parameter.
        internal_coordination_scopes = {
            ct.DATASET.value: dataset_scope
        }
        if coordination_scopes is not None:
            internal_coordination_scopes.update(coordination_scopes)
        vcv = VitessceConfigView(
            component_str, internal_coordination_scopes, x, y, w, h)

        # Use the mapping parameter if component is scatterplot and the mapping is not None
        if mapping is not None:
            [et_scope] = self.add_coordination(ct.EMBEDDING_TYPE)
            et_scope.set_value(mapping)
            vcv.use_coordination(et_scope)

        if isinstance(props, dict):
            vcv.set_props(**props)

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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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
            assert isinstance(c_type, ct) or isinstance(c_type, str)
            if isinstance(c_type, str):
                c_type_str = c_type
            else:
                c_type_str = c_type.value
            prev_scopes = list(self.config["coordinationSpace"][c_type_str].keys(
            )) if c_type_str in self.config["coordinationSpace"].keys() else []
            scope = VitessceConfigCoordinationScope(
                c_type_str, _get_next_scope(prev_scopes))
            if scope.c_type not in self.config["coordinationSpace"]:
                self.config["coordinationSpace"][scope.c_type] = {}
            self.config["coordinationSpace"][scope.c_type][scope.c_scope] = scope
            result.append(scope)
        return result

    def set_coordination_value(self, c_type, c_scope, c_value):
        """
        Set the value for a coordination scope. If a coordination object for the coordination type does not yet exist in the coordination space, it will be created.

        :param str c_type: The coordination type for this coordination scope.
        :param str c_scope: The coordination scope name.
        :param any c_value: The value for the coordination scope. Optional.

        :returns: The coordination scope instance.
        :rtype: VitessceConfigCoordinationScope
        """
        scope = VitessceConfigCoordinationScope(c_type, c_scope, c_value)
        if scope.c_type not in self.config["coordinationSpace"]:
            self.config["coordinationSpace"][scope.c_type] = {}
        self.config["coordinationSpace"][scope.c_type][scope.c_scope] = scope
        return scope

    def link_views(self, views, c_types, c_values=None):
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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            vc.layout(hconcat(v1, vconcat(v2, v3)))

        .. code-block:: python
            :emphasize-lines: 8

            from vitessce import VitessceConfig, Component as cm

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v2 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            v3 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
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
                        x_min + (w / num_views) * i,
                        x_min + (w / num_views) * (i + 1),
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
                        y_min + (h / num_views) * i,
                        y_min + (h / num_views) * (i + 1),
                    )

        # Recursively set the values (x,y,w,h) for each view.
        _layout(view_concat, 0, 12, 0, 12)

        # TODO: decide how to handle views that were omitted from the `view_concat` parameter
        # TODO: decide how to handle .add_view() being called after .layout() has been called

        return self

    def to_dict(self, base_url=None):
        """
        Convert the view config instance to a dict object.

        :param str base_url: Optional parameter for non-remote data to specify the url from which the data will be served.

        :returns: The view config as a dict. Useful for serializing to JSON.
        :rtype: dict
        """
        return {
            **self.config,
            "datasets": [d.to_dict(base_url) for d in self.config["datasets"]],
            "coordinationSpace": dict([
                (c_type, dict([
                    (c_scope_name, c_scope.c_value) for c_scope_name, c_scope in c_scopes.items()
                ])) for c_type, c_scopes in self.config["coordinationSpace"].items()
            ]),
            # TODO: compute the x,y,w,h values if not explicitly defined
            "layout": [c.to_dict() for c in self.config["layout"]]
        }

    def get_routes(self):
        """
        Convert the routes for this view config from the datasets.

        :returns: A list of server routes.
        :rtype: list[starlette.routing.Route]
        """
        routes = []
        for d in self.config["datasets"]:
            routes += d.get_routes()
        return routes

    def to_python(self):
        """
        Convert the VitessceConfig instance to a one-line Python code snippet that can be used to generate it.

        :returns: (A list of classes from the vitessce package used in the code block, The formatted code block)
        :rtype: (list[str], str)
        """
        classes_to_import = OrderedDict()
        classes_to_import[VitessceChainableConfig.__name__] = True
        code_block = f'{VitessceChainableConfig.__name__}({make_params_repr(self._to_py_params())})'

        for vcd in self.config["datasets"]:
            vcd_file_list_contents = ', '.join(
                [repr(f) for f in vcd._get_files()])
            vcd_obj_list_contents = ', '.join(
                [repr(f) for f in vcd._get_objects()])
            add_dataset_func = self.add_dataset.__name__
            add_dataset_params_list = [
                make_params_repr(vcd._to_py_params()),
            ]
            if len(vcd._get_files()) > 0:
                add_dataset_params_list.append(
                    f'files=[{vcd_file_list_contents}]')
                classes_to_import[VitessceConfigDatasetFile.__name__] = True
            if len(vcd._get_objects()) > 0:
                add_dataset_params_list.append(
                    f'objs=[{vcd_obj_list_contents}]')
            add_dataset_params = ', '.join(add_dataset_params_list)
            code_block += f'.{add_dataset_func}({add_dataset_params})'
            for obj in vcd._get_objects():
                if "vitessce" in sys.modules and obj.__class__.__name__ in dict(inspect.getmembers(sys.modules["vitessce"])):
                    classes_to_import[obj.__class__.__name__] = True
        for c_type, c_obj in self.config["coordinationSpace"].items():
            if c_type != ct.DATASET.value:
                for c_scope_name, c_scope in c_obj.items():
                    set_coordination_func = self.set_coordination_value.__name__
                    set_coordination_params = make_params_repr(
                        c_scope._to_py_params())
                    code_block += f'.{set_coordination_func}({set_coordination_params})'

        for vcv in self.config["layout"]:
            dataset_for_view = self.get_dataset_by_coordination_scope_name(
                vcv.get_coordination_scope(ct.DATASET.value))
            if dataset_for_view is not None:
                dataset_uid = dataset_for_view.get_uid()
            elif len(self.config["datasets"]) == 1:
                # If there is only one dataset available, assume it is the dataset for this view.
                dataset_uid = self.config["datasets"][0].get_uid()
            else:
                raise ValueError(
                    "At least one dataset must be present in the config before adding a view.")
            add_view_params_dict = {
                "dataset_uid": dataset_uid,
            }
            add_view_params_dict.update(vcv._to_py_params())
            if vcv.get_props() is not None:
                add_view_params_dict["props"] = vcv.get_props()
            add_view_func = self.add_view.__name__
            add_view_params = make_params_repr(add_view_params_dict)
            code_block += f'.{add_view_func}({add_view_params})'
        formatted_code_block = black.format_str(
            code_block, mode=black.FileMode())
        return list(classes_to_import), formatted_code_block

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
        vc = VitessceConfig(
            name=config["name"], description=config["description"], schema_version=config["version"])

        # Add each dataset from the incoming config.
        for d in config["datasets"]:
            new_dataset = vc.add_dataset(uid=d["uid"], name=d["name"])
            for f in d["files"]:
                new_dataset.add_file(
                    url=f.get("url"),
                    data_type=f["type"],
                    file_type=f["fileType"],
                    options=f.get("options")
                )
        if 'coordinationSpace' in config:
            for c_type in config['coordinationSpace'].keys():
                if c_type != ct.DATASET.value:
                    c_obj = config['coordinationSpace'][c_type]
                    vc.config['coordinationSpace'][c_type] = {}
                    for c_scope_name, c_scope_value in c_obj.items():
                        scope = VitessceConfigCoordinationScope(
                            c_type, c_scope_name)
                        scope.set_value(c_scope_value)
                        vc.config['coordinationSpace'][c_type][c_scope_name] = scope

        for c in config['layout']:
            c_coord_scopes = c['coordinationScopes'] if 'coordinationScopes' in c.keys() else {
            }
            if len(config["datasets"]) > 1 and ct.DATASET.value not in c_coord_scopes:
                raise ValueError(
                    "Multiple datasets are present, so every view must have an explicit dataset coordination scope.")
            new_view = VitessceConfigView(
                c['component'], c_coord_scopes, c['x'], c['y'], c['w'], c['h'])
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

        # The data object may modify the view config if it implements the auto_view_config() method.
        obj.auto_view_config(vc)

        return vc

    def widget(self, **kwargs):
        """
        Instantiate a VitessceWidget object based on this config.

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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            vc.layout(v1)
            vw = vc.widget()
            vw
        """
        from .widget import VitessceWidget  # TODO: Move import back to top when this is factored out.
        return VitessceWidget(self, **kwargs)

    def web_app(self, **kwargs):
        """
        Launch the http://vitessce.io web app using this config.

        :param str theme: The theme name, either "light" or "dark". By default, "auto", which selects light or dark based on operating system preferences.
        :param int port: The port to use when serving data objects on localhost. By default, 8000.
        :param base_url: If the web app is being accessed remotely (i.e. the data is being served from a remote machine), specify the base URL here. If serving and accessing the data on the same machine, keep as None to use a localhost URL.
        :type base_url: str or None
        :param bool open: Should the browser be opened to the web app URL? By default, True.

        :returns: The URL of the web app (containing the Vitessce configuration as URL-encoded JSON).
        :rtype: str

        .. code-block:: python
            :emphasize-lines: 6

            from vitessce import VitessceConfig, Component as cm, CoordinationType as ct

            vc = VitessceConfig()
            my_dataset = vc.add_dataset(name='My Dataset')
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            vc.layout(v1)
            vc.web_app()
        """
        from .widget import launch_vitessce_io  # TODO: Move import back to top when this is factored out.
        return launch_vitessce_io(self, **kwargs)

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
            v1 = vc.add_view(cm.SPATIAL, dataset=my_dataset)
            vc.layout(v1)

            config_dict = vc.export(to="S3")
        """
        from .export import (export_to_s3, export_to_files)  # TODO: Move import back to top when this is factored out.
        if to == "S3":
            return export_to_s3(self, *args, **kwargs)
        elif to == "files":
            return export_to_files(self, *args, **kwargs)
        else:
            raise ValueError("Unknown export destination.")


class VitessceChainableConfig(VitessceConfig):
    """
    A class to represent a Vitessce view config, where the methods ``add_dataset``, ``add_view``, and ``set_coordination_value`` return self (the config instance). This class inherits from ``VitessceConfig``.
    """

    def __init__(self, **kwargs):
        """
        Construct a Vitessce view config object.

        :param \*\*kwargs:  Takes the same arguments as the constructor on the ``VitessceConfig`` class.

        .. code-block:: python
            :emphasize-lines: 3

            from vitessce import VitessceChainableConfig

            vc = VitessceChainableConfig(name='My Config')
        """
        super().__init__(**kwargs)

    def __copy__(self):
        new_vc = VitessceChainableConfig(
            name=self.config["name"], description=self.config["description"], schema_version=self.config["version"])
        new_vc.config = self.config.copy()
        return new_vc

    def add_dataset(self, copy=True, **kwargs):
        """
        Add a dataset to this config.

        :param \*\*kwargs: Takes the same arguments as the ``add_dataset`` method on the ``VitessceConfig`` class.

        :returns: The config instance.
        :rtype: VitessceChainableConfig
        """
        if copy:
            new_vc = copy_module.copy(self)
            return new_vc.add_dataset(copy=False, **kwargs)
        super().add_dataset(**kwargs)
        return self

    def add_view(self, component, copy=True, **kwargs):
        """
        Add a view to this config.

        :param component: Takes the same arguments as the ``add_view`` method on the ``VitessceConfig`` class.
        :param \*\*kwargs: Takes the same arguments as the ``add_view`` method on the ``VitessceConfig`` class.

        :returns: The config instance.
        :rtype: VitessceChainableConfig
        """
        if copy:
            new_vc = copy_module.copy(self)
            return new_vc.add_view(component, copy=False, **kwargs)
        super().add_view(component, **kwargs)
        return self

    def set_coordination_value(self, c_type, c_scope, c_value, copy=True):
        """
        Add a coordination value to this config.

        :param c_type: Takes the same arguments as the ``set_coordination_value`` method on the ``VitessceConfig`` class.
        :param c_scope: Takes the same arguments as the ``set_coordination_value`` method on the ``VitessceConfig`` class.
        :param c_value: Takes the same arguments as the ``set_coordination_value`` method on the ``VitessceConfig`` class.

        :returns: The config instance.
        :rtype: VitessceChainableConfig
        """
        if copy:
            new_vc = copy_module.copy(self)
            return new_vc.set_coordination_value(c_type, c_scope, c_value, copy=False)
        super().set_coordination_value(c_type, c_scope, c_value)
        return self
