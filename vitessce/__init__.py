import sys
from warnings import warn

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
