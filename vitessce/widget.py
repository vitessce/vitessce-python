import importlib.util
from urllib.parse import quote_plus
import json

# Widget dependencies
import anywidget
from traitlets import Unicode, Dict, Int, Bool
import time

# Server dependencies
import asyncio
from hypercorn.config import Config
from hypercorn.asyncio import serve

from starlette.applications import Starlette
from starlette.middleware import Middleware
from starlette.middleware.cors import CORSMiddleware
from threading import Thread
import socket

MAX_PORT_TRIES = 1000
DEFAULT_PORT = 8000


def run_server_loop(app, port):
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    config = Config()
    config.bind = [f"localhost:{port}"]

    # As of Hypercorn 0.11.0, need to explicitly set signal handlers to a no-op
    # (otherwise it will try to set signal handlers assuming it is on the main thread which throws an error)
    loop.run_until_complete(
        serve(app, config, shutdown_trigger=lambda: asyncio.Future()))
    loop.close()


def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0


def get_base_url_and_port(port, next_port, proxy=False, base_url=None):
    if port is None:
        use_port = next_port
        next_port += 1
        port_tries = 1
        while is_port_in_use(use_port) and port_tries < MAX_PORT_TRIES:
            use_port = next_port
            next_port += 1
            port_tries += 1
    else:
        use_port = port

    if base_url is None:
        if proxy:
            if importlib.util.find_spec('jupyter_server_proxy') is None:
                raise ValueError(
                    "To use the widget through a proxy, jupyter-server-proxy must be installed.")
            base_url = f"proxy/{use_port}"
        else:
            base_url = f"http://localhost:{use_port}"

    return base_url, use_port, next_port


def serve_routes(routes, use_port):
    if len(routes) > 0:
        middleware = [
            Middleware(CORSMiddleware, allow_origins=[
                       '*'], allow_methods=["OPTIONS", "GET"], allow_headers=['Range'])
        ]
        app = Starlette(debug=True, routes=routes, middleware=middleware)

        t = Thread(target=run_server_loop, args=(app, use_port))
        t.start()
        time.sleep(1)


def launch_vitessce_io(config, theme='light', port=None, base_url=None, open=True):
    import webbrowser
    base_url, use_port, _ = get_base_url_and_port(
        port, DEFAULT_PORT, base_url=base_url)
    config_dict = config.to_dict(base_url=base_url)
    routes = config.get_routes()
    serve_routes(routes, use_port)
    vitessce_url = f"http://vitessce.io/?theme={theme}&url=data:," + quote_plus(
        json.dumps(config_dict))
    if open:
        webbrowser.open(vitessce_url)
    return vitessce_url


ESM = """
import * as d3 from "https://esm.sh/d3-require@1.3.0";
import React from 'https://unpkg.com/es-react@16.13.1/react.js';
import ReactDOM from 'https://unpkg.com/es-react@16.13.1/react-dom.js';

function asEsModule(component) {
  return {
    __esModule: true,
    default: component,
  };
}

const e = React.createElement;

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

export function render(view) {
    const jsPackageVersion = view.model.get('js_package_version');
    let customRequire = d3.require;
    const customJsUrl = view.model.get('custom_js_url');
    if(customJsUrl.length > 0) {
        customRequire = d3.requireFrom(async () => {
            return customJsUrl;
        });
    }

    const aliasedRequire = customRequire.alias({
        "react": React,
        "react-dom": ReactDOM
    });

    const Vitessce = React.lazy(() => aliasedRequire(`vitessce@${jsPackageVersion}`).then(vitessce => asEsModule(vitessce.Vitessce)));

    function VitessceWidget(props) {
        const { model } = props;

        const config = prependBaseUrl(model.get('config'), model.get('proxy'));
        const height = model.get('height');
        const theme = model.get('theme') === 'auto' ? (prefersDark ? 'dark' : 'light') : model.get('theme');

        const divRef = React.useRef();

        React.useEffect(() => {
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

        const onConfigChange = React.useCallback((config) => {
            model.set('config', config);
            model.save_changes();
        }, [model]);

        const vitessceProps = { height, theme, config, onConfigChange };

        return e('div', { ref: divRef, style: { height: height + 'px' } },
            e(React.Suspense, { fallback: e('div', {}, 'Loading...') },
                e(Vitessce, vitessceProps)
            )
        );
    }

    ReactDOM.render(e(VitessceWidget, { model: view.model }), view.el);
}
"""


class VitessceWidget(anywidget.AnyWidget):
    """
    A class to represent a Jupyter widget for Vitessce.
    """
    _module = Unicode(ESM).tag(sync=True)

    # Widget specific property.
    # Widget properties are defined as traitlets. Any property tagged with `sync=True`
    # is automatically synced to the frontend *any* time it changes in Python.
    # It is synced back to Python from the frontend *any* time the model is touched.
    config = Dict({}).tag(sync=True)
    height = Int(600).tag(sync=True)
    theme = Unicode('auto').tag(sync=True)
    proxy = Bool(False).tag(sync=True)

    next_port = DEFAULT_PORT

    js_package_version = Unicode('2.0.2').tag(sync=True)
    custom_js_url = Unicode('').tag(sync=True)

    def __init__(self, config, height=600, theme='auto', port=None, proxy=False, js_package_version='2.0.2', custom_js_url=''):
        """
        Construct a new Vitessce widget.

        :param config: A view config instance.
        :type config: VitessceConfig
        :param str theme: The theme name, either "light" or "dark". By default, "auto", which selects light or dark based on operating system preferences.
        :param int height: The height of the widget, in pixels. By default, 600.
        :param int port: The port to use when serving data objects on localhost. By default, 8000.
        :param bool proxy: Is this widget being served through a proxy, for example with a cloud notebook (e.g. Binder)?

        .. code-block:: python
            :emphasize-lines: 4

            from vitessce import VitessceConfig, VitessceWidget

            vc = VitessceConfig.from_object(my_scanpy_object)
            vw = vc.widget()
            vw
        """

        base_url, use_port, VitessceWidget.next_port = get_base_url_and_port(
            port, VitessceWidget.next_port, proxy=proxy)
        config_dict = config.to_dict(base_url=base_url)
        routes = config.get_routes()

        super(VitessceWidget, self).__init__(
            config=config_dict, height=height, theme=theme, proxy=proxy,
            js_package_version=js_package_version, custom_js_url=custom_js_url
        )

        serve_routes(routes, use_port)

    def _get_coordination_value(self, coordination_type, coordination_scope):
        obj = self.config['coordinationSpace'][coordination_type]
        obj_scopes = list(obj.keys())
        if coordination_scope is not None:
            if coordination_scope in obj_scopes:
                return obj[coordination_scope]
            else:
                raise ValueError(
                    f"The specified coordination scope '{coordination_scope}' could not be found for the coordination type '{coordination_type}'. Known coordination scopes are {obj_scopes}")
        else:
            if len(obj_scopes) == 1:
                auto_coordination_scope = obj_scopes[0]
                return obj[auto_coordination_scope]
            elif len(obj_scopes) > 1:
                raise ValueError(
                    f"The coordination scope could not be automatically determined because multiple coordination scopes exist for the coordination type '{coordination_type}'. Please specify one of {obj_scopes} using the scope parameter.")
            else:
                raise ValueError(
                    f"No coordination scopes were found for the coordination type '{coordination_type}'.")

    def get_cell_selection(self, scope=None):
        return self._get_coordination_value('cellSelection', scope)
