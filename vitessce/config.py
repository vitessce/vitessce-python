import json
from uuid import uuid4

from .constants import (
    CoordinationType as ct,
    Component as cm,
    DataType as dt,
    FileType as ft
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
    def __init__(self, uid, name, files):
        if files is not None:
            for f in files:
                assert type(f) == VitessceConfigDatasetFile
        
        self.dataset = {
            "uid": uid,
            "name": name or "",
            "files": files or [],
        }
    
    def add_file(self, url, data_type, file_type):
        """
        Add a new file definition to this dataset instance.

        Parameters
        ----------
        url : str

        data_type : str or DataType Enum

        file_type : str or FileType Enum

        Returns
        -------

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

    def to_dict(self):
        return {
            **self.dataset,
            "files": [ f.to_dict() for f in self.dataset["files"] ],
        }

class VitessceConfigView:
    """
    A class to represent a view (i.e. visualization component) in the Vitessce view config layout.
    """
    def __init__(self, component, coordination_scopes, x, y, w, h):
        self.view = {
            "component": component,
            "coordinationScopes": coordination_scopes,
            "x": 0,
            "y": 0,
            "w": 1,
            "h": 1
        }
    
    def use_coordination(self, *c_scopes):
        """
        Attach a coordination scope to this view instance. All views using the same coordination scope for a particular coordination type will effectively be linked together.

        Parameters
        ----------
        *c_scopes : VitessceConfigCoordinationScope
            Variable number of coordination scope instances can be passed.

        Returns
        -------
        VitessceConfigView
            The view instance, to allow chaining.
        """
        for c_scope in c_scopes:
            assert type(c_scope) == VitessceConfigCoordinationScope
            self.view["coordinationScopes"][c_scope.c_type] = c_scope.c_scope
        return self
    
    def to_dict(self):
        return self.view

class VitessceConfigCoordinationScope:
    """
    A class to represent a coordination scope in the Vitessce view config coordination space.
    """
    def __init__(self, c_type, c_scope):
        self.c_type = c_type
        self.c_scope = c_scope
        self.c_value = None
    
    def set_value(self, v):
        self.c_value = v
        return self

class VitessceConfig:
    """
    A class to represent a Vitessce view config.
    """

    def __init__(self, config=None, name=None, description=None):
        """
        Construct a Vitessce view config object.

        Parameters
        ----------
        config : dict
            Optionally provide an existing Vitessce view config as a dict to allow manipulation through the VitessceConfig API.
        name : str
            A name for the view config. Optional.
        description : str
            A description for the view config. Optional.
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

        if config is None:
            # No existing config was provided.
            if name is None:
                self.config["name"] = ""
            else:
                self.config["name"] = Name
            if description is None:
                self.config["description"] = ""
            else:
                self.config["description"] = description
        
        else:
            # An existing config was provided.

            # TODO: Validate the incoming config.
            
            self.config["name"] = config["name"]
            self.config["description"] = config["description"]

            # Add each dataset from the incoming config.
            for d in config["datasets"]:
                new_dataset = self.add_dataset(uid=d["uid"], name=d["name"])
                for f in d["files"]:
                    new_file = new_dataset.add_file(
                        url=f["url"],
                        data_type=f["type"],
                        file_type=f["fileType"]
                    )
                        
            # TODO: Add each coordination scope from the incoming config.

            # TODO: Add the components (layout) from the incoming config.


    def add_dataset(self, name="", uid=None, files=None):
        """
        Add a dataset to the config.

        Parameters
        ----------
        name : str
            A name for this dataset.
        uid : str
            A unique identifier for this dataset. Optional. If None, will be automatically generated.
        files : list of VitessceConfigDatasetFile
            A list of file objects. Optional. Files can also be added to the dataset via the .add_file method on the returned dataset instance.

        Returns
        -------
        VitessceConfigDataset
            The instance for the new dataset.
        """
        uid = uid if uid is not None else _get_next_scope([ d.dataset['uid'] for d in self.config["datasets"] ])
        assert type(uid) == str
        vcd = VitessceConfigDataset(uid, name, files)
        self.config["datasets"].append(vcd)
        [d_scope] = self.add_coordination(ct.DATASET)
        d_scope.set_value(uid)
        return vcd
    
    def add_view(self, dataset, component, x=None, y=None, w=None, h=None, mapping=None):
        """
        Add a view to the config.

        Parameters
        ----------
        dataset : VitessceConfigDataset
            A dataset instance to be used for the data visualized in this view.
        component : str or Component Enum
            A component name, either as a string or using the Component enum values.
        mapping : str
            An optional convenience parameter for setting the EMBEDDING_TYPE coordination scope value. Only applicable to the SCATTERPLOT component.
        x : int
            The horizontal position of the view. Must be an integer between 0 and 11. Optional.
        y : int
            The vertical position of the view. Must be an integer between 0 and 11. Optional.
        w : int
            The width of the view. Must be an integer between 1 and 12. Optional.
        h : int
            The height of the view. Must be an integer between 1 and 12. Optional.

        Returns
        -------
        VitessceConfigView
            The instance for the new view.
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
            for scope_name, dataset_scope in self.config["coordinationSpace"]["dataset"].items()
            if dataset_scope.c_value == dataset.dataset["uid"]
        ] if "dataset" in self.config["coordinationSpace"].keys() else []
        if len(dataset_matches) == 1:
            dataset_scope = dataset_matches[0]
        else:
            raise ValueError("The dataset parameter could not be found in the coordination space.")

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

        Parameters
        ----------
        *c_types : str or CoordinationType Enum
            Variable number of coordination types can be passed.
        
        Returns
        -------
        list of VitessceConfigCoordinationScope
            The instances for the new scope objects corresponding to each coordination type. These can be linked to views via the .use_coordination function on the VitessceConfigView class.
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
        
    def to_dict(self):
        """
        Convert the view config instance to a dict object.

        Returns
        -------
        dict
            The view config as a dict. Useful if serialization to JSON is required.
        """
        return {
            **self.config,
            "datasets": [ d.to_dict() for d in self.config["datasets"] ],
            "coordinationSpace": dict([
                (c_type, dict([
                    (c_scope_name, c_scope.c_value) for c_scope_name, c_scope in c_scopes.items() 
                ])) for c_type, c_scopes in self.config["coordinationSpace"].items()
            ]),
            "layout": [ c.to_dict() for c in self.config["layout"] ]
        }



