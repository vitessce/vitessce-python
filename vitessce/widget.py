import ipywidgets as widgets
from traitlets import Unicode, Dict

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
    config = Dict({
        "version": "0.1.0",
        "description": "High Bit Depth (uint16) Multiplex Immunofluorescence Imaging",
        "layers": [
        {
            "name": "raster",
            "type": "RASTER",
            "fileType": "raster.json",
            "url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/spraggins/spraggins.raster.json"
        }
        ],
        "name": "Spraggins",
        "public": True,
        "staticLayout": [
        {
            "component": "spatial",
            "props": {
            "view": {
                "zoom": -6.5,
                "target": [
                20000,
                20000,
                0
                ]
            }
            },
            "x": 0,
            "y": 0,
            "w": 9,
            "h": 2
        },
        {
            "component": "layerController",
            "x": 9,
            "y": 0,
            "w": 3,
            "h": 2
        }
        ]
    }).tag(sync=True)
