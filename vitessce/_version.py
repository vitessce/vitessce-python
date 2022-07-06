from pathlib import Path
from json import loads

# Module version
py_version_info = (1, 0, 9)

_js_info = loads((Path(__file__).parent.parent / 'js/package.json').read_text())
js_version_info = _js_info['version'].split('.')

# Module version accessible using vitessce.__version__
__version__ = '%s.%s.%s' % (
    py_version_info[0], py_version_info[1], py_version_info[2])
