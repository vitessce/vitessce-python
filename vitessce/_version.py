# The version of the Python package (`vitessce` on PyPI).
py_version_info = (1, 0, 9)
# The version of the JS package, should match the value in js/package.json (`vitessce-jupyter` on NPM).
js_version_info = (0, 1, 14)
# The maximum view config schema version supported by the bundled Vitessce JS package js/package.json.
vc_version_info = (1, 0, 7)

# Module version accessible using vitessce.__version__
__version__ = '%s.%s.%s' % (
    py_version_info[0], py_version_info[1], py_version_info[2])
