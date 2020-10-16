import json
from uuid import uuid4

from .constants import CoordinationTypes as ct

def get_next_scope(prev_scopes):
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

"""
function getNextScope(prevScopes) {
  // Keep an ordered list of valid characters.
  const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  // Store the value of the next character for each position
  // in the new string.
  // For example, [0] -> "A", [1] -> "B", [0, 1] -> "AB"
  const nextCharIndices = [0];

  // Generate a new scope name,
  // potentially conflicting with an existing name.
  // Reference: https://stackoverflow.com/a/12504061
  function next() {
    const r = [];
    nextCharIndices.forEach((charIndex) => {
      r.unshift(chars[charIndex]);
    });
    let increment = true;
    for (let i = 0; i < nextCharIndices.length; i++) {
      const val = ++nextCharIndices[i];
      if (val >= chars.length) {
        nextCharIndices[i] = 0;
      } else {
        increment = false;
        break;
      }
    }
    if (increment) {
      nextCharIndices.push(0);
    }
    return r.join('');
  }

  let nextScope;
  do {
    nextScope = next();
  } while (prevScopes.includes(nextScope));
  return nextScope;
}
"""

class VitessceConfigDatasetFile:
    def __init__(self, url, data_type, file_type):
        self.file = {
            "url": url,
            "type": data_type,
            "fileType": file_type,
        }
    
    def to_dict(self):
        return self.file

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

    def to_dict(self):
        return {
            **self.dataset,
            "files": [ f.to_dict() for f in self.dataset["files"] ],
        }

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
    def __init__(self, config=None, name=None, description=None):
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

        if config is None:
            if name is None:
                self.config["name"] = ""
            else:
                self.config["name"] = Name
            if description is None:
                self.config["description"] = ""
            else:
                self.config["description"] = description
        
        else:
            # TODO: validate the incoming config

            self.config["name"] = config["name"]
            self.config["description"] = config["description"]

            for d in config["datasets"]:
                new_dataset = self.add_dataset(uid=d["uid"], name=d["name"])
                for f in d["files"]:
                    new_file = new_dataset.add_file_url(
                        url=f["url"],
                        data_type=f["type"],
                        file_type=f["fileType"]
                    )


    def add_dataset(self, **kwargs):
        kwargs['uid'] = kwargs['uid'] if 'uid' in kwargs else get_next_scope([ d.dataset['uid'] for d in self.config["datasets"] ])
        vcd = VitessceConfigDataset(**kwargs)
        self.config["datasets"].append(vcd)
        self.config["coordinationSpace"]["dataset"][get_next_scope(list(self.config["coordinationSpace"]["dataset"].keys()))] = kwargs['uid']
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
            scope = VitessceConfigCoordinationScope(c_type.value, get_next_scope(list(self.config["coordinationSpace"][c_type.value].keys())))
            if scope.c_type not in self.config["coordinationSpace"]:
                self.config["coordinationSpace"][scope.c_type] = {}
            self.config["coordinationSpace"][scope.c_type][scope.c_scope] = scope
            result.append(scope)
        return result
        
    def to_dict(self):

        return {
            **self.config,
            "datasets": [ d.to_dict() for d in self.config["datasets"] ]
        }



