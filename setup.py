from __future__ import print_function
from setuptools import setup, find_packages
import os
from os.path import join as pjoin
from distutils import log

from jupyter_packaging import (
    create_cmdclass,
    install_npm,
    ensure_targets,
    combine_commands,
    get_version,
)


here = os.path.dirname(os.path.abspath(__file__))

log.set_verbosity(log.DEBUG)
log.info('setup.py entered')
log.info('$PATH=%s' % os.environ['PATH'])

name = 'vitessce'

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

# Get vitessce version
version = get_version(pjoin(name, '_version.py'))

js_dir = pjoin(here, 'js')

# Representative files that should exist after a successful build
jstargets = [
    pjoin(js_dir, 'dist', 'index.js'),
]

data_files_spec = [
    ('share/jupyter/nbextensions/vitessce-jupyter', 'vitessce/nbextension', '*.*'),
    ('share/jupyter/labextensions/vitessce-jupyter', 'vitessce/labextension', '**'),
    ('share/jupyter/labextensions/vitessce-jupyter', '.', 'install.json'),
    ('etc/jupyter/nbconfig/notebook.d', '.', 'vitessce-jupyter.json'),
]

cmdclass = create_cmdclass('jsdeps', data_files_spec=data_files_spec)
cmdclass['jsdeps'] = combine_commands(
    install_npm(js_dir, npm=['npm'],
                build_cmd='build'), ensure_targets(jstargets),
)


extras_require = {
    'building': [
        'build==0.1.0',
    ],
    'jupyter': [
        # Required for developing jupyter extensions.
        'jupyterlab==3.1.14',
        'jupyter_server==1.11.0',
        # Traitlets and markupsafe versions need to be pinned.
        # Reference: https://github.com/vitessce/vitessce-python/pull/143#discussion_r854340283
        'traitlets==4.3.3',
        'markupsafe==2.0.1',
    ],
    'testing': [
        'pytest>=6.2.4',
        'anndata==0.7.8',
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
        'ipywidgets<=7.7.2',
        'jupyterlab-widgets<=1.1.1',
        'hypercorn>=0.11.0',
        'ujson>=4.0.1',
        'starlette==0.14.0',

        # aiofiles is not explicitly referenced in our code,
        # but it is an implicit dependency of starlette==0.14.0.
        # https://github.com/encode/starlette/issues/49
        # Upgrading starlette will remove this dependency.
        'aiofiles>=0.6.0'
    ]
}

# Option for user to install all runtime deps.
extras_require['all'] = extras_require['proxy'] + extras_require['notebook']

# Option for developers to install all runtime deps + all development deps.
extras_require['dev'] = sum(extras_require.values(), [])

setup_args = dict(
    name=name,
    version=version,
    description='Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    include_package_data=True,
    install_requires=[
        'zarr>=2.5.0',
        'numcodecs>=0.5.7',
        'scipy>=1.2.1',
        'negspy>=0.2.24',
        'pandas>=1.1.2',
        'black>=21.11b1',
        'numpy>=1.21.2',
        'generate-tiff-offsets>=0.1.7',
    ],
    extras_require=extras_require,
    packages=find_packages(),
    zip_safe=False,
    cmdclass=cmdclass,
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
