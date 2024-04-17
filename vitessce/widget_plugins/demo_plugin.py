from esbuild_py import transform

plugin_esm = transform("""
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