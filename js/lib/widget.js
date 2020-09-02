import { DOMWidgetView, DOMWidgetModel } from '@jupyter-widgets/base';
import React from 'react';
import ReactDOM from 'react-dom';
import { Vitessce } from 'vitessce';
import 'vitessce/dist/es/production/static/css/index.css';
import './widget.css';

// See example.py for the kernel counterpart to this file.

const config = {
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
    "public": true,
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
};

class VitessceWidget extends React.Component {

    componentDidMount() {
        console.log("componentDidMount");
    }
    render() {
        return React.createElement('div', { },
            React.createElement(Vitessce, { config, height: 500, theme: 'dark' })
        );
    }
}

// Custom Model. Custom widgets models must at least provide default values
// for model attributes, including
//
//  - `_view_name`
//  - `_view_module`
//  - `_view_module_version`
//
//  - `_model_name`
//  - `_model_module`
//  - `_model_module_version`
//
//  when different from the base class.

// When serialiazing the entire widget state for embedding, only values that
// differ from the defaults will be specified.
export const VitessceModel = DOMWidgetModel.extend({
    defaults: _.extend(DOMWidgetModel.prototype.defaults(), {
        _model_name : 'VitessceModel',
        _view_name : 'VitessceView',
        _model_module : 'vitessce-jupyter',
        _view_module : 'vitessce-jupyter',
        _model_module_version : '0.1.0',
        _view_module_version : '0.1.0',
        value : 'Hello World!'
    })
});

export class VitessceView extends DOMWidgetView {
    render() {
        super.render();
        ReactDOM.render(
            this._render(this.model, this),
            this.el,
        );
    }

    _render(model, view) {
        return React.createElement(VitessceWidget, { model, view });
    }
}