from enum import Enum


BASE_URL_PLACEHOLDER = "{{ base_url }}"


# Reference: https://stackoverflow.com/a/50473952
class DocEnum(Enum):
    def __new__(cls, value, doc):
        self = object.__new__(cls)
        self._value_ = value
        self.__doc__ = doc
        return self


def norm_enum(enum_val, expected_enum_class=None):
    assert isinstance(enum_val, str) or isinstance(enum_val, expected_enum_class)
    # We don't actually use the expected_enum_class,
    # since it would not account for things like plugin coordination types, etc.
    # But we can pass it around anyway and in the future could use
    # it for extra checks, e.g. upon debug flags being true.
    if isinstance(enum_val, str):
        return enum_val
    else:
        return enum_val.value


class CoordinationType(DocEnum):
    """
    An enum type representing a coordination type in the Vitessce coordination model.
    The term coordination type refers to a parameter to be coordinated, and its programming-language-like type.
    For example, the ``SPATIAL_ZOOM`` coordination type represents a coordination of the zoom level of a spatial view, which can take a float value.
    """
    META_COORDINATION_SCOPES = "metaCoordinationScopes", "Shared representation of view-level coordinationScopes."
    META_COORDINATION_SCOPES_BY = "metaCoordinationScopesBy", "Shared representation of view-level coordinationScopesBy."
    DATASET = "dataset", "The identifier for the dataset associated with a view."
    OBS_TYPE = "obsType", "The type of entity represented by each observation."
    FEATURE_TYPE = "featureType", "The type of entity represented by each feature."
    FEATURE_VALUE_TYPE = "featureValueType", "The type of value stored in the observation-by-feature matrix."
    OBS_LABELS_TYPE = 'obsLabelsType', "Feature for displaying additional obs sets' data in heatmap/scatterplot/spatial tooltips."
    EMBEDDING_TYPE = "embeddingType", "The type of embedding used for a scatterplot view, such as PCA or t-SNE."
    EMBEDDING_ZOOM = "embeddingZoom", "The zoom level of an embedding scatterplot view."
    EMBEDDING_ROTATION = "embeddingRotation", "The rotation of an embedding scatterplot view."
    EMBEDDING_TARGET_X = "embeddingTargetX", "The x-coordinate of the center of an embedding scatterplot view."
    EMBEDDING_TARGET_Y = "embeddingTargetY", "The y-coordinate of the center of an embedding scatterplot view."
    EMBEDDING_TARGET_Z = "embeddingTargetZ", "The z-coordinate of the center of an embedding scatterplot view."
    EMBEDDING_OBS_SET_POLYGONS_VISIBLE = 'embeddingObsSetPolygonsVisible', "Whether polygon boundaries for each selected obsSet are visible in the embedding scatterplt."
    EMBEDDING_OBS_SET_LABELS_VISIBLE = 'embeddingObsSetLabelsVisible', "Whether labels for each selected obsSet are visible in the embedding scatterplot."
    EMBEDDING_OBS_SET_LABEL_SIZE = 'embeddingObsSetLabelSize', "The size of labels for selected obsSets in the embedding scatterplot."
    EMBEDDING_OBS_OPACITY = 'embeddingObsOpacity', "Opacity of cells in embedding (points in scatterplot view)"
    EMBEDDING_OBS_RADIUS = 'embeddingObsRadius', "Radius of cells in embedding (points in scatterplot view)"
    EMBEDDING_OBS_RADIUS_MODE = 'embeddingObsRadiusMode', "Radius mode of cells in embedding (points in scatterplot view) - auto or manual"
    EMBEDDING_OBS_OPACITY_MODE = 'embeddingObsOpacityMode', "Opacity mode of cells in embedding (points in scatterplot view) - auto or manual"
    EMBEDDING_CELL_OPACITY = 'embeddingCellOpacity', "Deprecated"
    EMBEDDING_CELL_RADIUS = 'embeddingCellRadius', "Deprecated"
    EMBEDDING_CELL_RADIUS_MODE = 'embeddingCellRadiusMode', "Deprecated"
    EMBEDDING_CELL_OPACITY_MODE = 'embeddingCellOpacityMode', "Deprecated"
    SPATIAL_ZOOM = "spatialZoom", "The zoom level of a spatial view."
    SPATIAL_ROTATION_X = 'spatialRotationX', "The x rotation of a 3d spatial view."
    SPATIAL_ROTATION_Y = 'spatialRotationY', "The y rotation of a 3d spatial view."
    SPATIAL_ROTATION_Z = 'spatialRotationZ', "The z rotation of a 3d spatial view."
    SPATIAL_ROTATION_ORBIT = 'spatialRotationOrbit', "The rotation orbit in degrees of a 3d spatial view."
    SPATIAL_ORBIT_AXIS = 'spatialOrbitAxis', "The orbital axis of a 3d spatial view."
    SPATIAL_AXIS_FIXED = 'spatialAxisFixed', "Boolean for whether or not the target axis of a spatial view is fixed."
    SPATIAL_TARGET_X = "spatialTargetX", "The x-coordinate of the center of a spatial view."
    SPATIAL_TARGET_Y = "spatialTargetY", "The y-coordinate of the center of a spatial view."
    SPATIAL_TARGET_Z = "spatialTargetZ", "The z-coordinate of the center of a spatial view."
    HEATMAP_ZOOM_X = "heatmapZoomX", "The x-axis zoom level of a heatmap view."
    HEATMAP_ZOOM_Y = "heatmapZoomY", "The y-axis zoom level of a heatmap view."
    HEATMAP_TARGET_X = "heatmapTargetX", "The x-coordinate of the center of a heatmap view."
    HEATMAP_TARGET_Y = "heatmapTargetY", "The y-coordinate of the center of a heatmap view."
    OBS_FILTER = "obsFilter", "A subset of cells to include after filtering."
    OBS_HIGHLIGHT = "obsHighlight", "A subset of cells to highlight."
    OBS_SET_SELECTION = "obsSetSelection", "A subset of cell sets to select."
    OBS_SET_HIGHLIGHT = "obsSetHighlight", "A subset of cell sets to highlight."
    OBS_SET_COLOR = "obsSetColor", "A mapping from cell sets to colors."
    CELL_FILTER = "cellFilter", "Deprecated"
    CELL_HIGHLIGHT = "cellHighlight", "Deprecated"
    CELL_SET_SELECTION = "cellSetSelection", "Deprecated"
    CELL_SET_HIGHLIGHT = "cellSetHighlight", "Deprecated"
    CELL_SET_COLOR = "cellSetColor", "Deprecated"
    FEATURE_FILTER = "featureFilter", "A subset of genes to include after filtering."
    FEATURE_HIGHLIGHT = "featureHighlight", "A subset of genes to highlight."
    FEATURE_SELECTION = "featureSelection", "A subset of genes to select."
    FEATURE_VALUE_COLORMAP = "featureValueColormap", "The colormap to use for the gene expression scale."
    FEATURE_VALUE_TRANSFORM = 'featureValueTransform', "Function to use to transform feature values."
    FEATURE_VALUE_COLORMAP_RANGE = "featureValueColormapRange", "The range of gene expression values to map."
    GENE_FILTER = "geneFilter", "Deprecated"
    GENE_HIGHLIGHT = "geneHighlight", "Deprecated"
    GENE_SELECTION = "geneSelection", "Deprecated"
    GENE_EXPRESSION_COLORMAP = "geneExpressionColormap", "Deprecated"
    GENE_EXPRESSION_COLORMAP_RANGE = "geneExpressionColormapRange", "Deprecated"
    OBS_COLOR_ENCODING = "obsColorEncoding", "The color encoding to use for cell entities."
    CELL_COLOR_ENCODING = "cellColorEncoding", "Deprecated"
    SPATIAL_IMAGE_LAYER = 'spatialImageLayer', "Layer definitions for the imagery in the spatial view."
    SPATIAL_SEGMENTATION_LAYER = 'spatialSegmentationLayer', "Layer definitions for the segmentations in the spatial view."
    SPATIAL_POINT_LAYER = 'spatialPointLayer', "Layer definitions for the points in the spatial view."
    SPATIAL_NEIGHBORHOOD_LAYER = 'spatialNeighborhoodLayer', "Layer definitions for the neighborhoods in the spatial view."
    SPATIAL_RASTER_LAYERS = 'spatialRasterLayers', "Deprecated"
    SPATIAL_CELLS_LAYER = 'spatialCellsLayer', "Deprecated"
    SPATIAL_MOLECULES_LAYER = 'spatialMoleculesLayer', "Deprecated"
    SPATIAL_NEIGHBORHOODS_LAYER = 'spatialNeighborhoodsLayer', "Deprecated"
    GENOMIC_ZOOM_X = "genomicZoomX", "The zoom level of a higlass view, X dimension."
    GENOMIC_ZOOM_Y = "genomicZoomY", "The zoom level of a higlass view, Y dimension."
    GENOMIC_TARGET_X = "genomicTargetX", "The x-coordinate of the center of a higlass view."
    GENOMIC_TARGET_Y = "genomicTargetY", "The y-coordinate of the center of a higlass view."
    ADDITIONAL_CELL_SETS = "additionalCellSets", "Deprecated"
    ADDITIONAL_OBS_SETS = "additionalObsSets", "User-defined cell sets."
    MOLECULE_HIGHLIGHT = 'moleculeHighlight', "Deprecated"
    GATING_FEATURE_SELECTION_X = 'gatingFeatureSelectionX', "Feature for the x-axis of the gating scatterplot."
    GATING_FEATURE_SELECTION_Y = 'gatingFeatureSelectionY', "Feature for the y-axis of the gating scatterplot."
    FEATURE_VALUE_TRANSFORM_COEFFICIENT = 'featureValueTransformCoefficient', "Coefficient to transform values in the gating scatterplot."
    TOOLTIPS_VISIBLE = 'tooltipsVisible', "Boolean for whether or not tooltips are visible, used by the scatterplot, spatial, and heatmap views."


class ViewType(DocEnum):
    """
    An enum type representing a view type in the visualization layout.
    """
    SCATTERPLOT = "scatterplot", "The scatterplot component can be used for visualization of 2-dimensional embeddings."
    SPATIAL = "spatial", "The spatial component can be used for visualization of cells, molecules, or images in spatial coordinates."
    DESCRIPTION = "description", "The description component can display short informational text about a dataset."
    STATUS = "status", "The status component can display contextual information such as hover states or error messages."
    HEATMAP = "heatmap", "The heatmap component can be used to view a cell by gene expression matrix."
    LAYER_CONTROLLER = "layerController", "The layer controller can be used to manipulate channel settings of the images rendered by the spatial component."
    GENOMIC_PROFILES = "genomicProfiles", "The higlass component can be used to visualize genome-wide ATAC-seq profiles."
    OBS_SETS = "obsSets", "Observation sets"
    OBS_SET_SIZES = "obsSetSizes", "Observation set sizes bar plot"
    OBS_SET_FEATURE_VALUE_DISTRIBUTION = "obsSetFeatureValueDistribution", "Violin plot visualizing the distribution of feature values per observation set."
    FEATURE_VALUE_HISTOGRAM = "featureValueHistogram", "Histogram visualizing the distribution of values for a selected feature across all observations."
    FEATURE_LIST = "featureList", "The feature list selector"
    GATING = "gating", "A gating scatterplot"
    CELL_SETS = "cellSets", "Deprecated"
    CELL_SET_SIZES = "cellSetSizes", "Deprecated"
    CELL_SET_EXPRESSION = "cellSetExpression", "Deprecated"
    GENES = "genes", "Deprecated"


class DataType(DocEnum):
    """
    An enum type representing the type of data contained in a file.
    """
    OBS_LABELS = "obsLabels", "Alternate label for each observation."
    OBS_EMBEDDING = "obsEmbedding", "Embedding coordinates for each observation."
    OBS_LOCATIONS = "obsLocations", "Spatial coordinates for each observation."
    OBS_FEATURE_MATRIX = "obsFeatureMatrix", "Matrix of feature (column) values for each observation (row)."
    OBS_SETS = "obsSets", "Sets of observations, representing clusters or classes."
    FEATURE_LABELS = "featureLabels", "Alternate label for each feature."
    IMAGE = "image", "Image"
    OBS_SEGMENTATIONS = "obsSegmentations", "Segmentation for each observation."
    NEIGHBORHOODS = "neighborhoods", "The spatial cell neighborhoods data type."
    GENOMIC_PROFILES = "genomic-profiles", "The genomic profiles data type, used by HiGlass 1D quantitative tracks."
    CELLS = "cells", "Deprecated"
    CELL_SETS = "cell-sets", "Deprecated"
    EXPRESSION_MATRIX = "expression-matrix", "Deprecated"
    MOLECULES = "molecules", "Deprecated"
    RASTER = "raster", "Deprecated"


class FileType(DocEnum):
    """
    An enum type representing the file format or schema to which a file conforms.
    """
    ANNDATA_ZARR = "anndata.zarr", "Joint file type for AnnData objects"
    ANNDATA_H5AD = "anndata.h5ad", "Joint file type for AnnData objects"
    OBS_EMBEDDING_CSV = 'obsEmbedding.csv', "File type for obsEmbedding values stored in a CSV file"
    OBS_LOCATIONS_CSV = 'obsLocations.csv', "File type for obsLocations values stored in a CSV file"
    OBS_LABELS_CSV = 'obsLabels.csv', "File type for obsLabels values stored in a CSV file"
    FEATURE_LABELS_CSV = 'featureLabels.csv', "File type for featureLabels values stored in a CSV file"
    OBS_FEATURE_MATRIX_CSV = 'obsFeatureMatrix.csv', "File type for obsFeatureMatrix stored in a CSV file"
    OBS_SEGMENTATIONS_JSON = 'obsSegmentations.json', "File type for obsSegmentations polygons stored in a JSON file"
    OBS_SETS_CSV = 'obsSets.csv', "File type for obsSets stored in a CSV file"
    OBS_SETS_JSON = 'obsSets.json', "File type for obsSets stored in a JSON file"
    IMAGE_OME_ZARR = "image.ome-zarr", "File type for images stored as OME-NGFF Zarr stores."
    OBS_FEATURE_MATRIX_ANNDATA_ZARR = 'obsFeatureMatrix.anndata.zarr', "File type for obsFeatureMatrix stored in an AnnData object saved to a Zarr store"
    OBS_SETS_ANNDATA_ZARR = 'obsSets.anndata.zarr', "File type for obsSets stored in an AnnData object saved to a Zarr store"
    OBS_EMBEDDING_ANNDATA_ZARR = 'obsEmbedding.anndata.zarr', "File type for obsEmbedding values stored in an AnnData object saved to a Zarr store"
    OBS_LOCATIONS_ANNDATA_ZARR = 'obsLocations.anndata.zarr', "File type for obsLocations values stored in an AnnData object saved to a Zarr store"
    OBS_SEGMENTATIONS_ANNDATA_ZARR = 'obsSegmentations.anndata.zarr', "File type for obsSegmentations polygons stored in an AnnData object saved to a Zarr store"
    OBS_LABELS_ANNDATA_ZARR = 'obsLabels.anndata.zarr', "File type for obsLabels stored in an AnnData object saved to a Zarr store"
    FEATURE_LABELS_ANNDATA_ZARR = 'featureLabels.anndata.zarr', "File type for featureLabels stored in an AnnData object saved to a Zarr store"
    EXPRESSION_MATRIX_ZARR = "expression-matrix.zarr", "The Zarr-based expression matrix file type."
    CELLS_JSON = "cells.json", "The JSON-based cells file type."
    MOLECULES_JSON = "molecules.json", "The JSON-based molecules file type."
    NEIGHBORHOODS_JSON = "neighborhoods.json", "The JSON-based neighborhoods file type."
    RASTER_JSON = "raster.json", "The JSON-based raster manifest file type."
    CELL_SETS_JSON = "cell-sets.json", "The JSON-based cell sets file type."
    CLUSTERS_JSON = "clusters.json", "A JSON-based expression matrix file type (this file type is poorly named)."
    GENES_JSON = "genes.json", "A JSON-based expression matrix file type."
    GENOMIC_PROFILES_ZARR = "genomic-profiles.zarr", "The Zarr-based genomic profile (multivec) file type."
    ANNDATA_CELLS_ZARR = "anndata-cells.zarr", "The Zarr-based cells file type from an anndata object."
    ANNDATA_CELL_SETS_ZARR = "anndata-cell-sets.zarr", "The Zarr-based cell-sets file type from an anndata object."
    ANNDATA_EXPRESSION_MATRIX_ZARR = "anndata-expression-matrix.zarr", "The Zarr-based expression matrix file type from an anndata object."
    OBS_SEGMENTATIONS_CELLS_JSON = "obsSegmentations.cells.json", "The JSON-based cells file type for obsSegmentations."
    OBS_LOCATIONS_CELLS_JSON = "obsLocations.cells.json", "The JSON-based cells file type for obsLocations."
    OBS_EMBEDDING_CELLS_JSON = "obsEmbedding.cells.json", "The JSON-based cells file type for obsEmbedding."
    OBS_SETS_CELL_SETS_JSON = "obsSets.cell-sets.json", "The JSON-based cell sets file type for obsSets."
