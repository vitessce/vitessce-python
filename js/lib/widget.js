import { DOMWidgetView, DOMWidgetModel } from '@jupyter-widgets/base';
import React, { useEffect, useRef, useCallback } from 'react';
import ReactDOM from 'react-dom';
import extend from 'lodash/extend';
import { Vitessce } from 'vitessce';
import packageJson from '../package.json';
import 'vitessce/dist/es/production/static/css/index.css';
import './widget.css';

// See widget.py for the kernel counterpart to this file.

const prefersDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;

function VitessceWidget(props) {
  const {
    model
  } = props;

  const config = model.get('config');
  const height = model.get('height');
  const theme = model.get('theme') === 'auto' ? (prefersDark ? 'dark' : 'light') : model.get('theme');

  const divRef = useRef();

  useEffect(() => {
    if(!divRef.current) {
      return () => {};
    }

    function handleMouseEnter() {
      const jpn = divRef.current.closest('.jp-Notebook');
      if(jpn) {
        jpn.style.overflow = "hidden";
      }
    }
    function handleMouseLeave() {
      const jpn = divRef.current.closest('.jp-Notebook');
      if(jpn) {
        jpn.style.overflow = "auto";
      }
    }
    divRef.current.addEventListener("mouseenter", handleMouseEnter);
    divRef.current.addEventListener("mouseleave", handleMouseLeave);

    return () => {
      divRef.current.removeEventListener("mouseenter", handleMouseEnter);
      divRef.current.removeEventListener("mouseleave", handleMouseLeave);
    };
  }, [divRef]);

  const onConfigChange = useCallback((config) => {
    model.set('config', config);
    model.save_changes();
  }, [model]);

  return React.createElement('div', { className: 'vitessce-widget', ref: divRef, style: { height: `${height}px` } },
    React.createElement(Vitessce, { config, onConfigChange, height, theme }),
  );
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
    defaults: extend(DOMWidgetModel.prototype.defaults(), {
        _model_name : 'VitessceModel',
        _view_name : 'VitessceView',
        _model_module : 'vitessce-jupyter',
        _view_module : 'vitessce-jupyter',
        _model_module_version : packageJson.version,
        _view_module_version : packageJson.version,
        config : {},
        height: 600,
        theme: 'auto',
    }),
});

export class VitessceView extends DOMWidgetView {
    render() {
        super.render();
        ReactDOM.render(
            this._render(this.model, this),
            this.el,
        );
        
        setTimeout(() => {
          window.dispatchEvent(new Event('resize'));
        }, 500);
    }

    _render(model, view) {
        return React.createElement(VitessceWidget, { model, view });
    }

    remove() {
      ReactDOM.unmountComponentAtNode(this.el);
      return super.remove();
    }
}