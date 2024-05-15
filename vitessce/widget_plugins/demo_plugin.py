from esbuild_py import transform

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
    function DemoView(props) {
        const { coordinationScopes } = props;
        const [{
            obsType,
        }, {
            setObsType,
        }] = useCoordination(['obsType'], coordinationScopes);

        return (
            <div className="demo">
                <p>Demo plugin view</p>
                <p>obsType: {obsType}</p>
                <button>This is an example button</button>
            </div>
        );
    }

    const pluginViewTypes = [
        new PluginViewType('demo', DemoView, ['obsType']),
    ];
    return { pluginViewTypes };
}
export default { createPlugins };
""")
"""
Example of a minimal plugin view that gets the obsType coordination value from the coordination space and renders a button.
This plugin view is not meant to be useful for end-users, but rather to demonstrate how to develop a plugin view that uses coordination (and uses eslint_py for JSX transformation).

:meta hide-value:

.. code-block:: python

    from vitessce.widget_plugins import demo_plugin_esm

    # ...
    vc.widget(plugin_esm=demo_plugin_esm)
"""
