# Widget dependencies
import ipywidgets as widgets
from traitlets import Unicode, Dict, Int
import time

# Server dependencies
import asyncio
from hypercorn.config import Config
from hypercorn.asyncio import serve

from starlette.applications import Starlette
from starlette.middleware import Middleware
from starlette.middleware.cors import CORSMiddleware
from threading import Thread

# Config creation dependencies
from .routes import create_obj_routes
from .config import VitessceConfig

# See js/lib/widget.js for the frontend counterpart to this file.

def run_server_loop(app, port):
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    config = Config()
    config.bind = [f"localhost:{port}"]

    # As of Hypercorn 0.11.0, need to explicitly set signal handlers to a no-op
    # (otherwise it will try to set signal handlers assuming it is on the main thread which throws an error)
    loop.run_until_complete(serve(app, config, shutdown_trigger=lambda: asyncio.Future()))
    loop.close()
    

@widgets.register
class VitessceWidget(widgets.DOMWidget):
    """
    A class to represent a Jupyter widget for Vitessce.
    """

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

    next_port = 8000

    def __init__(self, config, height=600, theme='dark', port=None):
        """
        Construct a new Vitessce widget.

        :param config: A view config instance.
        :type config: VitessceConfig
        :param str theme: The theme name, either "light" or "dark". By default, "dark".
        :param int height: The height of the widget, in pixels. By default, 600.
        :param int port: The port to use when serving data objects on localhost. By default, 8000.

        .. code-block:: python
            :emphasize-lines: 4

            from vitessce import VitessceConfig, VitessceWidget

            vc = VitessceConfig.from_object(my_scanpy_object)
            vw = VitessceWidget(vc)
            vw
        """

        if port is None:
            use_port = VitessceWidget.next_port
            VitessceWidget.next_port += 1
        else:
            use_port = port

        assert type(config) == VitessceConfig
        
        routes = []
        def on_obj(obj, dataset_uid, obj_i):
            obj_file_defs, obj_routes = create_obj_routes(obj, use_port, dataset_uid, obj_i)
            for obj_route in obj_routes:
                routes.append(obj_route)
            return obj_file_defs
        config_dict = config.to_dict(on_obj=on_obj)

        super(VitessceWidget, self).__init__(config=config_dict, height=height, theme=theme)
        
        if len(routes) > 0:
            middleware = [
                Middleware(CORSMiddleware, allow_origins=['*'])
            ]
            app = Starlette(debug=True, routes=routes, middleware=middleware)
            
            t = Thread(target=run_server_loop, args=(app, use_port))
            t.start()
            time.sleep(1)
            
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

