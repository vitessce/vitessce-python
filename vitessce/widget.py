# Widget dependencies
import ipywidgets as widgets
from traitlets import Unicode, Dict, Int

# Server dependencies
import asyncio
from hypercorn.config import Config
from hypercorn.asyncio import serve

from starlette.applications import Starlette

# Config creation dependencies
from .config import create_config_and_routes

# See js/lib/widget.js for the frontend counterpart to this file.

@widgets.register
class VitessceWidget(widgets.DOMWidget):
    """An example widget."""

    # Name of the widget view class in front-end
    _view_name = Unicode('VitessceView').tag(sync=True)

    # Name of the widget model class in front-end
    _model_name = Unicode('VitessceModel').tag(sync=True)

    # Name of the front-end module containing widget view
    _view_module = Unicode('vitessce-jupyter').tag(sync=True)

    # Name of the front-end module containing widget model
    _model_module = Unicode('vitessce-jupyter').tag(sync=True)

    # Version of the front-end module containing widget view
    _view_module_version = Unicode('^0.1.0').tag(sync=True)
    # Version of the front-end module containing widget model
    _model_module_version = Unicode('^0.1.0').tag(sync=True)
    
    # Widget specific property.
    # Widget properties are defined as traitlets. Any property tagged with `sync=True`
    # is automatically synced to the frontend *any* time it changes in Python.
    # It is synced back to Python from the frontend *any* time the model is touched.
    config = Dict({}).tag(sync=True)
    height = Int(600).tag(sync=True)
    theme = Unicode('dark').tag(sync=True)
    
    def __init__(self, **kwargs):
        
        routes = []
        if 'config' not in kwargs and 'data' in kwargs:
            kwargs['config'], routes = create_config_and_routes(kwargs['data'])

        super(VitessceWidget, self).__init__(**kwargs)
        
        if len(routes) > 0:
            app = Starlette(debug=True, routes=routes)
            
            # We cannot use asyncio.run() directly
            # since Jupyter runs in its own asyncio loop.
            # Reference: https://stackoverflow.com/a/61331974
            try:
                loop = asyncio.get_running_loop()
            except RuntimeError:
                loop = None
            
            if loop and loop.is_running():
                # As expected, there is already an event loop running:
                # the Jupyter event loop.
                task = loop.create_task(serve(app, Config()))
                task.add_done_callback(lambda t: print("Task done"))
            else:
                print('Error: did not find the expected notebook asyncio loop')
    
    def _get_coordination_value(self, coordination_type, coordination_scope):
        obj = self.config['coordinationSpace'][coordination_type]
        obj_scopes = list(obj.keys())
        if coordination_scope != None:
            if coordination_scope in obj_scopes:
                return obj[coordination_scope]
            else:
                raise ValueError(f"The specified coordination scope '{coordination_scope}' could not be found for the coordination type '{coordination_type}'. Known coordination scopes are {obj_scopes}")
        else:
            if len(obj_scopes) == 1:
                auto_coordination_scope = obj_scopes[0]
                return obj[auto_coordination_scope]
            elif len(obj_scopes) > 1:
                raise ValueError(f"The coordination scope could not be automatically determined because multiple coordination scopes exist for the coordination type '{coordination_type}'. Please specify one of {obj_scopes} using the scope parameter.")
            else:
                raise ValueError(f"No coordination scopes were found for the coordination type '{coordination_type}'.")
    
    def get_cell_selection(self, scope=None):
        return self._get_coordination_value('cellSelection', scope)

