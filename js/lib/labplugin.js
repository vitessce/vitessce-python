var plugin = require('./index');
var base = require('@jupyter-widgets/base');

module.exports = {
  id: 'vitessce-jupyter',
  requires: [base.IJupyterWidgetRegistry],
  activate: function(app, widgets) {
      widgets.registerWidget({
          name: 'vitessce-jupyter',
          version: plugin.version,
          exports: plugin
      });
  },
  autoStart: true
};

