import { DOMWidgetView, DOMWidgetModel } from '@jupyter-widgets/base';
import React, { useEffect, useRef, useCallback } from 'react';
import ReactDOM from 'react-dom';
import extend from 'lodash/extend';
import { Vitessce } from 'vitessce';
import packageJson from '../package.json';
import 'vitessce/dist/esm/index.css';
import './widget.css';

// See widget.py for the kernel counterpart to this file.

const prefersDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;

// The jupyter server may be running through a proxy,
// which means that the client needs to prepend the part of the URL before /proxy/8000 such as
// https://hub.gke2.mybinder.org/user/vitessce-vitessce-python-swi31vcv/proxy/8000/A/0/cells
function prependBaseUrl(config, proxy) {
  if(!proxy) {
    return config;
  }
  const { origin } = new URL(window.location.href);
  let baseUrl;
  const jupyterLabConfigEl = document.getElementById('jupyter-config-data');

  if (jupyterLabConfigEl) {
    // This is jupyter lab
    baseUrl = JSON.parse(jupyterLabConfigEl.textContent || '').baseUrl;
  } else {
    // This is jupyter notebook
    baseUrl = document.getElementsByTagName('body')[0].getAttribute('data-base-url');
  }
  return {
    ...config,
    datasets: config.datasets.map(d => ({
      ...d,
      files: d.files.map(f => ({
        ...f,
        url: `${origin}${baseUrl}${f.url}`,
      })),
    })),
  };
}

function VitessceWidget(props) {
  const {
    model
  } = props;

  const config = prependBaseUrl(model.get('config'), model.get('proxy'));
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
    function handleMouseLeave(event) {
      if(event.relatedTarget === null || (event.relatedTarget && event.relatedTarget.closest('.jp-Notebook')?.length)) return;
      const jpn = divRef.current.closest('.jp-Notebook');
      if(jpn) {
        jpn.style.overflow = "auto";
      }
    }
    divRef.current.addEventListener("mouseenter", handleMouseEnter);
    divRef.current.addEventListener("mouseleave", handleMouseLeave);

    return () => {
      if(divRef.current) {
        divRef.current.removeEventListener("mouseenter", handleMouseEnter);
        divRef.current.removeEventListener("mouseleave", handleMouseLeave);
      }
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
        proxy: false,
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
