# Widget dependencies
import ipywidgets as widgets
from traitlets import Unicode, Dict, Int

# Server dependencies
import asyncio
from hypercorn.config import Config
from hypercorn.asyncio import serve

from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route

# See js/lib/widget.js for the frontend counterpart to this file.

async def homepage(req):
    return JSONResponse({'hello': 'world'})

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
    height = Int().tag(sync=True)
    theme = Unicode().tag(sync=True)
    
    def __init__(self, **kwargs):
        super(VitessceWidget, self).__init__(**kwargs)
        
        config = kwargs['config']
        
        app = Starlette(debug=True, routes=[
            Route('/', homepage)
        ])
        
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

