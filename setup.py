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
LONG_DESCRIPTION = 'Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce'

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
    install_npm(js_dir, npm=['npm'], build_cmd='build'), ensure_targets(jstargets),
)

setup_args = dict(
    name=name,
    version=version,
    description='Jupyter widget facilitating interactive visualization of spatial single-cell data with Vitessce',
    long_description=LONG_DESCRIPTION,
    include_package_data=True,
    install_requires=[
        'ipywidgets>=7.6.0',
        'hypercorn>=0.11.0',
        'ujson>=4.0.1',
        'aiofiles>=0.6.0',
        'starlette==0.14.0',
        'zarr>=2.5.0',
        'numcodecs>=0.5.7',
        'scipy>=1.2.1',
        'negspy>=0.2.24',
        'generate-tiff-offsets>=0.1.7',
        'pandas>=1.1.2'
    ],
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
    ],
)

setup(**setup_args)
