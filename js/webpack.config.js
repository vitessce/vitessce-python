const path = require('path');
const fs = require('fs');
const version = require('./package.json').version;
const PnpWebpackPlugin = require('pnp-webpack-plugin');

const mode = 'development';
process.env.NODE_ENV = mode;

// Make sure any symlinks in the project folder are resolved:
function resolveApp(relativePath) {
    const appDirectory = fs.realpathSync(process.cwd());
    return path.resolve(appDirectory, relativePath);
}

// Custom webpack rules are generally the same for all webpack bundles, hence
// stored in a separate local variable.
const rules = [
    {
        test: /\.js$/,
        exclude: /node_modules/,
        use: {
            loader: require.resolve('babel-loader'),
            options: {
                customize: require.resolve(
                  'babel-preset-react-app/webpack-overrides'
                ),
            }
        },
    },
    {
        test: /\.css$/,
        use: [
            {
                loader: require.resolve('style-loader'),
            },
            {
                loader: require.resolve('css-loader'),
            },
        ],
    },
];

const resolve = {
    alias: {
      'txml/txml': 'txml/dist/txml'
    },
};

const resolveLoader = {
    plugins: [
      // Also related to Plug'n'Play, but this time it tells webpack to load its loaders
      // from the current package.
      PnpWebpackPlugin.moduleLoader(module),
    ],
};

module.exports = [
    {
        // Notebook extension
        //
        // This bundle only contains the part of the JavaScript that is run on
        // load of the notebook. This section generally only performs
        // some configuration for requirejs, and provides the legacy
        // "load_ipython_extension" function which is required for any notebook
        // extension.
        //
        mode,
        entry: './lib/extension.js',
        output: {
            filename: 'extension.js',
            path: path.resolve(__dirname, '..', 'vitessce', 'nbextension'),
            libraryTarget: 'amd'
        },
        devtool: 'cheap-source-map',
        module: {
            rules: rules
        },
        resolve,
        resolveLoader,
    },
    {
        // Bundle for the notebook containing the custom widget views and models
        //
        // This bundle contains the implementation for the custom widget views and
        // custom widget.
        // It must be an amd module
        //
        mode,
        entry: './lib/index.js',
        output: {
            filename: 'index.js',
            path: path.resolve(__dirname, '..', 'vitessce', 'nbextension'),
            libraryTarget: 'amd'
        },
        devtool: 'cheap-source-map',
        module: {
            rules: rules
        },
        externals: ['@jupyter-widgets/base'],
        resolve,
        resolveLoader,
    },
    {
        // Embeddable vitessce-jupyter bundle
        //
        // This bundle is generally almost identical to the notebook bundle
        // containing the custom widget views and models.
        //
        // The only difference is in the configuration of the webpack public path
        // for the static assets.
        //
        // It will be automatically distributed by unpkg to work with the static
        // widget embedder.
        //
        // The target bundle is always `dist/index.js`, which is the path required
        // by the custom widget embedder.
        //
        mode,
        entry: './lib/embed.js',
        output: {
            filename: 'index.js',
            path: path.resolve(__dirname, 'dist'),
            libraryTarget: 'amd',
            publicPath: 'https://unpkg.com/vitessce-jupyter@' + version + '/dist/'
        },
        devtool: 'cheap-source-map',
        module: {
            rules: rules
        },
        externals: ['@jupyter-widgets/base'],
        resolve,
        resolveLoader,
    },
    {
        // Jupyter Lab plugin
        mode,
        entry: './lib/labplugin.js',
        output: {
            filename: 'labplugin.js',
            path: path.resolve(__dirname, 'dist'),
            libraryTarget: 'amd'
        },
        devtool: 'cheap-source-map',
        module: {
            rules: rules
        },
        externals: ['@jupyter-widgets/base'],
        resolve,
        resolveLoader,
    },
];
