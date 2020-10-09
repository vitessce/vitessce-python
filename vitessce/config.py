import json
from uuid import uuid4

from constants import CoordinationTypes as ct


def next_uuid():
    return str(uuid4())

class VitessceConfigDatasetFile:
    def __init__(self, url, data_type, file_type):
        self.file = {
            "url": url,
            "type": data_type,
            "fileType": file_type,
        }

class VitessceConfigDataset:
    def __init__(self, uid=None, name="", files=[]):
        self.dataset = {
            "uid": uid,
            "name": name,
            "files": []
        }
    
    def add_file_url(self, url, **kwargs):
        self.dataset["files"].append(VitessceConfigDatasetFile(url=url, **kwargs))
        return self


class VitessceConfigView:
    def __init__(self, component, coordination_scopes):
        self.view = {
            "component": component,
            "coordinationScopes": coordination_scopes,
            "x": 0,
            "y": 0,
            "w": 1,
            "h": 1
        }
    
    def use_coordination(self, *c_scopes):
        for c_scope in c_scopes:
            assert type(c_scope) == VitessceConfigCoordinationScope
            self.view["coordinationScopes"][c_scope.c_type] = c_scope.c_scope
        return self

class VitessceConfigCoordinationScope:
    def __init__(self, c_type, c_scope):
        self.c_type = c_type
        self.c_scope = c_scope
        self.c_value = None

class VitessceConfig:
    def __init__(self, config=None, name="", description=""):
        # TODO: Allow passing an existing config, either as JSON or as a Python dict
        self.config = {
            "version": "1.0.0",
            "name": name,
            "description": description,
            "datasets": [],
            "coordinationSpace": {
                "dataset": {},
                "embeddingType": {},
            },
            "layout": [],
            "initStrategy": "auto"
        }

    def add_dataset(self, **kwargs):
        kwargs['uid'] = kwargs['uid'] if 'uid' in kwargs else next_uuid()
        vcd = VitessceConfigDataset(**kwargs)
        self.config["datasets"].append(vcd)
        self.config["coordinationSpace"]["dataset"][next_uuid()] = kwargs['uid']
        return vcd
    
    def add_view(self, dataset, component, mapping=None):
        assert type(dataset) == VitessceConfigDataset

        # Find the coordination scope name associated with the dataset
        dataset_matches = [
            scope
            for scope, dataset_id in self.config["coordinationSpace"]["dataset"].items()
            if dataset_id == dataset.dataset["uid"]
        ]
        if len(dataset_matches) == 1:
            dataset_scope = dataset_matches[0][0]
        else:
            raise ValueError("The dataset parameter could not be found in the coordination space.")

        coordination_scopes = {
            "dataset": dataset_scope,
        }
        vcv = VitessceConfigView(component=component, coordination_scopes=coordination_scopes)
        self.config["layout"].append(vcv)
        return vcv
    

    def add_coordination(self, *c_types):
        result = []
        for c_type in c_types:
            assert type(c_type) == ct
            scope = VitessceConfigCoordinationScope(c_type.value, next_uuid())
            if scope.c_type not in self.config["coordinationSpace"]:
                self.config["coordinationSpace"][scope.c_type] = {}
            self.config["coordinationSpace"][scope.c_type][scope.c_scope] = scope
            result.append(scope)
        return result
        
    def to_json(self):
        return self.config



