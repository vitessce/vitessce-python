import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Vitessce Widget Tutorial
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Vitessce custom plugin definition
        """
    )
    return


@app.cell
def _():
    from vitessce import (
        VitessceConfig,
        Component as cm,
        CoordinationType as ct,
        OmeTiffWrapper,
        MultiImageWrapper,
        VitesscePlugin
    )
    from oxc_py import transform
    return (
        MultiImageWrapper,
        OmeTiffWrapper,
        VitessceConfig,
        VitesscePlugin,
        cm,
        transform,
    )


@app.cell
def _(VitesscePlugin, transform):
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
    
        const CSS = `
            .chat {
                overflow-y: scroll;
            }
        `;
    
        function ChatView(props) {
        
            const [nextMessage, setNextMessage] = React.useState('');
            const [isLoading, setIsLoading] = React.useState(false);
            const [chatHistory, setChatHistory] = React.useState([]); // chatHistory is an array of message objects like [{ user, text }, ...]
                       
            async function handleClick() {     
                setChatHistory(prev => ([
                    ...prev,
                    { user: 'You', text: nextMessage },
                ]));
                setIsLoading(true);
                const [chatReceiveValue, chatReceiveBuffers] = await invokeCommand("chat_send", nextMessage, []);
                setChatHistory(prev => ([
                    ...prev,
                    { user: 'AI', text: chatReceiveValue.text },
                ]));
                setIsLoading(false);
            }
        
            return (
                <>
                    <style>{CSS}</style>
                    <div className="chat">
                        <p>Chat view</p>
                        <div>
                            {chatHistory.map(message => (
                                <p key={`${message.user}-${message.text}`}>
                                    <b>{message.user}</b>:
                                    <span>{message.text}</span>
                                </p>
                            ))}
                        </div>
                        <input type="text" value={nextMessage} onChange={e => setNextMessage(e.target.value)} disabled={isLoading} />
                        <button onClick={handleClick}>Send message</button>
                    </div>
                </>
            );
        }

        const pluginViewTypes = [
            new PluginViewType('chat', ChatView, []),
        ];
        return { pluginViewTypes };
    }
    export default { createPlugins };
    """)


    def handle_chat_message(message, buffers):
        return { "text": message.upper() }, []


    class ChatPlugin(VitesscePlugin):
        plugin_esm = PLUGIN_ESM
        commands = {
            "chat_send": handle_chat_message,
        }
    return (ChatPlugin,)


@app.cell
def _(MultiImageWrapper, OmeTiffWrapper, VitessceConfig, cm):
    vc = VitessceConfig(schema_version="1.0.15", name='Spraggins Multi-Modal', description='PAS + IMS + AF From https://portal.hubmapconsortium.org/browse/collection/6a6efd0c1a2681dc7d2faab8e4ab0bca')
    dataset = vc.add_dataset(name='Spraggins').add_object(
        MultiImageWrapper(
            image_wrappers=[
                OmeTiffWrapper(img_url='https://assets.hubmapconsortium.org/f4188a148e4c759092d19369d310883b/ometiff-pyramids/processedMicroscopy/VAN0006-LK-2-85-PAS_images/VAN0006-LK-2-85-PAS_registered.ome.tif?token=', name='PAS'),
                OmeTiffWrapper(img_url='https://assets.hubmapconsortium.org/2130d5f91ce61d7157a42c0497b06de8/ometiff-pyramids/processedMicroscopy/VAN0006-LK-2-85-AF_preIMS_images/VAN0006-LK-2-85-AF_preIMS_registered.ome.tif?token=', name='AF'),
                OmeTiffWrapper(img_url='https://assets.hubmapconsortium.org/be503a021ed910c0918842e318e6efa2/ometiff-pyramids/ometiffs/VAN0006-LK-2-85-IMS_PosMode_multilayer.ome.tif?token=', name='IMS Pos Mode'),
                OmeTiffWrapper(img_url='https://assets.hubmapconsortium.org/ca886a630b2038997a4cfbbf4abfd283/ometiff-pyramids/ometiffs/VAN0006-LK-2-85-IMS_NegMode_multilayer.ome.tif?token=', name='IMS Neg Mode')
            ],
            use_physical_size_scaling=True,
     )
    )
    spatial = vc.add_view(cm.SPATIAL, dataset=dataset)
    status = vc.add_view("chat", dataset=dataset)
    lc = vc.add_view(cm.LAYER_CONTROLLER, dataset=dataset).set_props(disableChannelsIfRgbDetected=True)
    vc.layout(spatial | (lc / status));
    return (vc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## 2. Create the Vitessce widget
        """
    )
    return


@app.cell
def _(ChatPlugin, vc):
    vw = vc.widget(plugins=[ChatPlugin()])
    vw
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
