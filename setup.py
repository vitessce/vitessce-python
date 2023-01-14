from __future__ import print_function
from setuptools import setup, find_packages
import os
from distutils import log

# Module version
py_version_info = (3, 0, 0)
__version__ = '%s.%s.%s' % (py_version_info[0], py_version_info[1], py_version_info[2])

# Setup
here = os.path.dirname(os.path.abspath(__file__))

log.set_verbosity(log.DEBUG)
log.info('setup.py entered')
log.info('$PATH=%s' % os.environ['PATH'])

name = 'vitessce'

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

extras_require = {
    'building': [
        'build==0.1.0',
    ],
    'testing': [
        'pytest>=6.2.4',
        'loompy>=3.0.6',
        'coverage>=6.3.2'
    ],
    'linting': [
        'flake8==3.8.4',
    ],
    'docs': [
        'sphinx==4.2.0',
        'sphinx-rtd-theme==1.0.0',
        'nbsphinx==0.8.8',
        'nbclean==0.3.2',
        # nbconvert and jinja2 versions need to be pinned.
        # Reference: https://github.com/vitessce/vitessce-python/issues/152
        'nbconvert==5.6.1',
        'jinja2==3.0.3',
    ],
    'proxy': [
        'jupyter-server-proxy>=1.5.2'
    ],
    'notebook': [
        # Needed only for notebook use:
        'anywidget==0.0.3',
        'uvicorn>=0.17.0',
        'ujson>=4.0.1',
        'starlette==0.14.0',
        'generate-tiff-offsets>=0.1.7',

        # aiofiles is not explicitly referenced in our code,
        # but it is an implicit dependency of starlette==0.14.0.
        # https://github.com/encode/starlette/issues/49
        # Upgrading starlette will remove this dependency.
        'aiofiles>=0.6.0'
    ],
}

# Option for user to install all runtime deps.
extras_require['all'] = extras_require['proxy'] + extras_require['notebook']

# Option for developers to install all runtime deps + all development deps.
extras_require['dev'] = sum(extras_require.values(), [])

setup_args = dict(
    name=name,
    version=__version__,
    description='Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    install_requires=[
        'zarr>=2.5.0',
        'numcodecs>=0.5.7',
        'scipy>=1.2.1',
        'negspy>=0.2.24',
        'pandas>=1.1.2',
        'black>=21.11b1',
        'numpy>=1.21.2',
        'anndata>=0.7.8,<0.9',
        'ome-zarr==0.2.1',
        'tifffile>=2020.10.1',
    ],
    extras_require=extras_require,
    packages=find_packages(),
    author='Gehlenborg Lab',
    author_email='',
    url='https://github.com/vitessce/vitessce-python',
    keywords=[
        'ipython',
        'jupyter',
        'widgets',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Framework :: IPython',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Multimedia :: Graphics',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)

setup(**setup_args)
