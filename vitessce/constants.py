from enum import Enum

# Reference: https://stackoverflow.com/a/50473952
class DocEnum(Enum):
    def __new__(cls, value, doc=None):
        self = object.__new__(cls)
        self._value_ = value
        if doc is not None:
            self.__doc__ = doc
        return self

class CoordinationType(DocEnum):
    """
    An enum type representing a coordination type in the Vitessce coordination model.
    The term coordination type refers to a parameter to be coordinated, and its programming-language-like type.
    For example, the ``SPATIAL_ZOOM`` coordination type represents a coordination of the zoom level of a spatial view, which can take a float value.
    """
    DATASET = "dataset", "The identifier for the dataset associated with a view."
    EMBEDDING_TYPE = "embeddingType", "The type of embedding used for a scatterplot view, such as PCA or t-SNE."
    EMBEDDING_ZOOM = "embeddingZoom", "The zoom level of an embedding scatterplot view."
    EMBEDDING_ROTATION = "embeddingRotation", "The rotation of an embedding scatterplot view."
    EMBEDDING_TARGET_X = "embeddingTargetX", "The x-coordinate of the center of an embedding scatterplot view."
    EMBEDDING_TARGET_Y = "embeddingTargetY", "The y-coordinate of the center of an embedding scatterplot view."
    EMBEDDING_TARGET_Z = "embeddingTargetZ", "The z-coordinate of the center of an embedding scatterplot view."
    SPATIAL_ZOOM = "spatialZoom"
    SPATIAL_ROTATION = "spatialRotation"
    SPATIAL_TARGET_X = "spatialTargetX"
    SPATIAL_TARGET_Y = "spatialTargetY"
    SPATIAL_TARGET_Z = "spatialTargetZ"
    HEATMAP_ZOOM_X = "heatmapZoomX"
    HEATMAP_ZOOM_Y = "heatmapZoomY"
    HEATMAP_TARGET_X = "heatmapTargetX"
    HEATMAP_TARGET_Y = "heatmapTargetY"
    CELL_FILTER = "cellFilter"
    CELL_SELECTION = "cellSelection"
    CELL_HIGHLIGHT = "cellHighlight"
    GENE_FILTER = "geneFilter"
    GENE_SELECTION = "geneSelection"
    GENE_HIGHLIGHT = "geneHighlight"
    GENE_EXPRESSION_COLORMAP = "geneExpressionColormap"
    GENE_EXPRESSION_COLORMAP_RANGE = "geneExpressionColormapRange"
    CELL_COLOR_ENCODING = "cellColorEncoding"
    SPATIAL_LAYERS = "spatialLayers"
    GENOMIC_ZOOM = "genomicZoom"
    GENOMIC_TARGET_X = "genomicTargetX"
    GENOMIC_TARGET_Y = "genomicTargetY"

class Component(DocEnum):
    SCATTERPLOT = "scatterplot"
    SPATIAL = "spatial"
    DESCRIPTION = "description"
    STATUS = "status"
    CELL_SETS = "cellSets"
    HEATMAP = "heatmap"
    LAYER_CONTROLLER = "layerController"
    HIGLASS = "higlass"
    CELL_SET_SIZES = "cellSetSizes"

class DataType(DocEnum):
    CELLS = "cells"
    CELL_SETS = "cell-sets"
    EXPRESSION_MATRIX = "expression-matrix"
    MOLECULES = "molecules"
    NEIGHBORHOODS = "neighborhoods"
    RASTER = "raster"

class FileType(DocEnum):
    EXPRESSION_MATRIX_ZARR = "expression-matrix.zarr"
    CLUSTERS_JSON = "clusters.json"
    GENES_JSON = "genes.json"
    CELLS_JSON = "cells.json"
    MOLECULES_JSON = "molecules.json"
    NEIGHBORHOODS_JSON = "neighborhoods.json"
    RASTER_JSON = "raster.json"
    CELL_SETS_JSON = "cell-sets.json"
