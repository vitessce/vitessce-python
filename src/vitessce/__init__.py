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

from .utils import (
    get_initial_coordination_scope_prefix,
    get_initial_coordination_scope_name,
    make_ids_csv_data_url,
    make_colors_csv_data_url,
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

from .wrappers import AbstractWrapper

# We allow installation without all of the dependencies that the widget requires.
# The imports below will fail in that case, and corresponding globals will be undefined.
try:
    from .widget import VitessceWidget, VitesscePlugin, data_server
except ModuleNotFoundError as e:  # pragma: no cover
    warn(f'Extra installs are necessary to use widgets: {e}')

try:
    from .wrappers import (
        OmeTiffWrapper,
        OmeZarrWrapper,
        MultiImageWrapper,
        CsvWrapper,
        JsonWrapper,
        AnnDataWrapper,
        MultivecZarrWrapper,
        ImageOmeTiffWrapper,
        ObsSegmentationsOmeTiffWrapper,
        ImageOmeZarrWrapper,
        ObsSegmentationsOmeZarrWrapper,
        SpatialDataWrapper,
        ObsSegmentationsNgPrecomputedWrapper,
        ObsPointsNgAnnotationsWrapper,
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
