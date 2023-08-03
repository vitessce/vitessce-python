import sys
from warnings import warn

from .config import (
    VitessceConfig,
    VitessceChainableConfig,
    VitessceConfigDatasetFile,
    hconcat,
    vconcat,
    CoordinationLevel,
)

from .repr import make_repr

from .constants import (
    CoordinationType,
    ViewType,
    DataType,
    FileType,
    # For backwards compatibility, also export ViewType as Component
    ViewType as Component,
    BASE_URL_PLACEHOLDER,
)

from .config_converter import (
    CellBrowserToAnndataZarrConverter,  # only exported for testing.
    convert_cell_browser_project_to_anndata,
)

from .wrappers import AbstractWrapper

# We allow installation without all of the dependencies that the widget requires.
# The imports below will fail in that case, and corresponding globals will be undefined.
try:
    from .widget import VitessceWidget, data_server
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use widgets: {e}')

try:
    from .wrappers import (
        OmeTiffWrapper,
        OmeZarrWrapper,
        MultiImageWrapper,
        CsvWrapper,
        AnnDataWrapper,
        MultivecZarrWrapper,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use wrappers: {e}')

try:
    from .export import (
        export_to_s3,
        export_to_files,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use exports: {e}')
