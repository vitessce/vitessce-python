from vitessce import (
    __version__, _version
)


def test_py_str_version():
    assert len(__version__.split('.')) == 3


def test_py_version():
    assert len(_version.py_version_info) == 3


def test_js_version():
    assert len(_version.js_version_info) == 3
