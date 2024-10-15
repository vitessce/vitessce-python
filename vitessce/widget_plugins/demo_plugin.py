from esbuild_py import transform
from ..widget import VitesscePlugin


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
        invokeCommand,
    } = utilsForPlugins;
    function DemoView(props) {
        const { coordinationScopes } = props;
        const [{
            obsType,
        }, {
            setObsType,
        }] = useCoordination(['obsType'], coordinationScopes);

        function handleClick() {
            console.log(invokeCommand('demo_command', "Hello from command", []));
        }

        return (
            <div className="demo">
                <p>Demo plugin view</p>
                <p>obsType: {obsType}</p>
                <button onClick={handleClick}>This is an example button</button>
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


def handle_demo_command(message, buffers):
    return message.upper(), []


class DemoPlugin(VitesscePlugin):
    """
    Example of a minimal plugin view that gets the obsType coordination value from the coordination space and renders a button.
    This plugin view is not meant to be useful for end-users, but rather to demonstrate how to develop a plugin view that uses coordination (and uses eslint_py for JSX transformation).

    .. code-block:: python

        from vitessce.widget_plugins import DemoPlugin

        # ...
        vc.widget(plugins=[DemoPlugin()])
    """
    plugin_esm = PLUGIN_ESM
    commands = {
        "demo_command": handle_demo_command,
    }
