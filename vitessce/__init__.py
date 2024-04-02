import sys
from warnings import warn

from .config import (
    CoordinationLevel,
    VitessceChainableConfig,
    VitessceConfig,
    VitessceConfigDatasetFile,
    hconcat,
    vconcat,
)
from .config_converter import (
    CellBrowserToAnndataZarrConverter,  # only exported for testing.
    convert_cell_browser_project_to_anndata,
)
from .constants import (
    BASE_URL_PLACEHOLDER,
    CoordinationType,
    DataType,
    FileType,
    ViewType,
)
from .constants import (
    # For backwards compatibility, also export ViewType as Component
    ViewType as Component,
)
from .repr import make_repr
from .utils import (
    get_initial_coordination_scope_name,
    get_initial_coordination_scope_prefix,
)
from .wrappers import AbstractWrapper

# We allow installation without all of the dependencies that the widget requires.
# The imports below will fail in that case, and corresponding globals will be undefined.
try:
    from .widget import VitessceWidget, data_server
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f"Extra installs are necessary to use widgets: {e}")

try:
    from .wrappers import (
        AnnDataWrapper,
        CsvWrapper,
        ImageOmeTiffWrapper,
        ImageOmeZarrWrapper,
        MultiImageWrapper,
        MultivecZarrWrapper,
        ObsSegmentationsOmeTiffWrapper,
        ObsSegmentationsOmeZarrWrapper,
        OmeTiffWrapper,
        OmeZarrWrapper,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f"Extra installs are necessary to use wrappers: {e}")

try:
    from .export import (
        export_to_files,
        export_to_s3,
    )
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f"Extra installs are necessary to use exports: {e}")
