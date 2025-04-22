Widget plugins
##############

Vitessce supports multiple types of `plugins <http://vitessce.io/docs/dev-plugins/>`_ defined using JavaScript code.

Leveraging concepts from `anywidget <https://github.com/manzt/anywidget>`_, we can define such plugins directly from Python: plugin developers can supply custom JavaScript code via a Python string.

The most minimal example of such plugin JavaScript code is the following:

.. code-block:: python

    PLUGIN_ESM = """
    function createPlugins(utilsForPlugins) {
        const {
            React,
            PluginFileType,
            PluginViewType,
            PluginCoordinationType,
            PluginJointFileType,
            z,
            useCoordination,
        } = utilsForPlugins;
        return {
            pluginViewTypes: undefined,
            pluginFileTypes: undefined,
            pluginCoordinationTypes: undefined,
            pluginJointFileTypes: undefined,
        };
    }
    export default { createPlugins };
    """


The plugin string must be defined as an EcmaScript Module (ESM) that exports a function named ``createPlugins``.
The ``createPlugins`` function is called (on the initial render of the Jupyter widget) with the ``utilsForPlugins`` argument (to facilitate dependency injection) and returns an object with the following properties:

- ``pluginViewTypes``: an array of objects that define the view types of the plugin.
- ``pluginFileTypes```: an array of objects that define the file types of the plugin.
- ``pluginCoordinationTypes``: an array of objects that define the coordination types of the plugin.
- ``pluginJointFileTypes``: an array of objects that define the joint file types of the plugin.

If defined, these plugin arrays are passed to the Vitessce component as `props <http://vitessce.io/docs/dev-plugins/>`_ with the same names.

**Note**: For maximum stability of plugins, we recommend that plugin developers document which version(s) of the vitessce Python package that plugins have been developed under.

--------------------------------
Passing plugin ESM to the widget
--------------------------------

The plugin string can be passed to the widget using the ``plugins`` parameter and passing a subclass of ``VitesscePlugin``:


.. code-block:: python

    from vitessce import VitessceConfig, VitesscePlugin

    class MyPlugin(VitesscePlugin):
        plugin_esm = PLUGIN_ESM

    vc = VitessceConfig(description="A Vitessce widget with a custom plugin")
    # Some more configuration here...

    plugin = MyPlugin()
    vc.widget(plugins=[plugin])


-------------------------------
Defining plugin views using JSX
-------------------------------

Vitessce plugin view types are defined as React components.
During typical React component development, JSX syntax is used.
However, JSX is not valid JavaScript and therefore must be transformed to valid JavaScript before it can be passed to the widget where it will be interpreted as ESM.
Vitessce plugin developers then have two options for defining React components for plugin view types:

* Use ``React.createElement`` directly (without JSX).
* Use the ``transform`` function from `oxc_py <https://github.com/keller-mark/oxc-py>`_ to perform JSX to JS transformation.

.. code-block:: python

    from oxc_py import transform

    PLUGIN_ESM = transform("""
    function createPlugins(utilsForPlugins) {
        const {
            React,
            PluginFileType,
            PluginViewType,
            PluginCoordinationType,
            PluginJointFileType,
            z,
            useCoordination,
        } = utilsForPlugins;

        function MyPluginView(props) {
            return (
                <p>Hello world from JSX!</p>
            );
        }

        const pluginViewTypes = [
            new PluginViewType('myPlugin', MyPluginView, []),
        ];
        return { pluginViewTypes };
    }
    export default { createPlugins };
    """)


To import additional dependencies, JavaScript (more specifically, ESM) can be dynamically imported from CDN (such as ``unpkg`` or ``esm.sh``) within the ``createPlugins`` function:

.. code-block:: python

    PLUGIN_ESM = """
    async function createPlugins(utilsForPlugins) {
        const {
            React,
            PluginFileType,
            PluginViewType,
            PluginCoordinationType,
            PluginJointFileType,
            z,
            useCoordination,
        } = utilsForPlugins;

        const d3 = await import('https://cdn.jsdelivr.net/npm/d3@7/+esm');

        // Do something with d3 here...

        return {
            pluginViewTypes: undefined,
            pluginFileTypes: undefined,
            pluginCoordinationTypes: undefined,
            pluginJointFileTypes: undefined,
        };
    }
    export default { createPlugins };
    """


To support more complex import scenarios, see `dynamic-importmap <https://github.com/keller-mark/dynamic-importmap>`_.



vitessce.widget_plugins
***********************

.. automodule:: vitessce.widget_plugins.demo_plugin
 :members:

.. automodule:: vitessce.widget_plugins.spatial_query
 :members: