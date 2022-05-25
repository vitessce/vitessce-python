import sys
from warnings import warn

from ._version import __version__

from .config import (
    VitessceConfig,
    VitessceChainableConfig,
    VitessceConfigDatasetFile,
    hconcat,
    vconcat,
)

from .repr import make_repr

from .constants import CoordinationType, Component, DataType, FileType

from .wrappers import AbstractWrapper


# We allow installation without all of the dependencies that the widget requires.
# The imports below will fail in that case, and corresponding globals will be undefined.
try:
    from .widget import VitessceWidget
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use widgets: {e}')

try:
    from .wrappers import (
        OmeTiffWrapper,
        MultiImageWrapper,
        AnnDataWrapper,
        SnapWrapper,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use wrappers: {e}')

try:
    from .entities import (
        CellSets,
        Cells,
        Molecules,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use entities: {e}')

try:
    from .export import (
        export_to_s3,
        export_to_files,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use exports: {e}')

try:
    if "google.colab" in sys.modules:  # pragma: no cover
        from google.colab import output

        output.enable_custom_widget_manager()
except ImportError:  # pragma: no cover
    pass


def _jupyter_labextension_paths():  # pragma: no cover
    """Called by Jupyter Lab Server to detect if it is a valid labextension and
    to install the widget

    Returns
    =======
    src: Source directory name to copy files from. Webpack outputs generated files
        into this directory and Jupyter Lab copies from this directory during
        widget installation
    dest: Destination directory name to install widget files to. Jupyter Lab copies
        from `src` directory into <jupyter path>/labextensions/<dest> directory
        during widget installation
    """
    return [{
        'src': 'labextension',
        'dest': 'vitessce-jupyter',
    }]


def _jupyter_nbextension_paths():  # pragma: no cover
    """Called by Jupyter Notebook Server to detect if it is a valid nbextension and
    to install the widget

    Returns
    =======
    section: The section of the Jupyter Notebook Server to change.
        Must be 'notebook' for widget extensions
    src: Source directory name to copy files from. Webpack outputs generated files
        into this directory and Jupyter Notebook copies from this directory during
        widget installation
    dest: Destination directory name to install widget files to. Jupyter Notebook copies
        from `src` directory into <jupyter path>/nbextensions/<dest> directory
        during widget installation
    require: Path to importable AMD Javascript module inside the
        <jupyter path>/nbextensions/<dest> directory
    """
    return [{
        'section': 'notebook',
        'src': 'nbextension',
        'dest': 'vitessce-jupyter',
        'require': 'vitessce-jupyter/extension'
    }]
