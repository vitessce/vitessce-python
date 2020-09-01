import { DOMWidgetView, JupyterLuminoWidget, DOMWidgetModel } from '@jupyter-widgets/base';
import React from 'react';
import ReactDOM from 'react-dom';

// See example.py for the kernel counterpart to this file.

class VitessceWidget extends React.Component {
    constructor(props) {
        super(props);
        this.el = React.createRef();
    }

    async componentDidMount() {
        const { model, view } = this.props;
        const childView = await view.create_child_view(model);
        /* JupyterPhosphorWidget.attach() expects this.el.current to be attached to the DOM. This
         * requires the view to be attached to the DOM. This is the case after view.displayed is
         * resolved.
         */
        await view.displayed;
        JupyterLuminoWidget.attach(childView.pWidget, this.el.current);
    }

    render() {
        return React.createElement('div', { ref: this.el }, 'Hello world!');
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