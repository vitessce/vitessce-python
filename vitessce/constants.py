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
    """
    An enum type representing a view type in the visualization layout.
    """
    SCATTERPLOT = "scatterplot", "The scatterplot component can be used for visualization of 2-dimensional embeddings."
    SPATIAL = "spatial", "The spatial component can be used for visualization of cells, molecules, or images in spatial coordinates."
    DESCRIPTION = "description", "The description component can display short informational text about a dataset."
    STATUS = "status", "The status component can display contextual information such as hover states or error messages."
    CELL_SETS = "cellSets", "The cell sets component can be used to view and manipulate hierarchical (or flat) sets of cells, including sets representing cell type clusters."
    HEATMAP = "heatmap", "The heatmap component can be used to view a cell by gene expression matrix."
    LAYER_CONTROLLER = "layerController", "The layer controller can be used to manipulate channel settings of the images rendered by the spatial component."
    HIGLASS = "higlass", "The higlass component can be used to visualize genome browser tracks and genome-wide interaction heatmaps."
    CELL_SET_SIZES = "cellSetSizes", "The cell set sizes component can display the quantities of cells in selected cell sets."

class DataType(DocEnum):
    """
    An enum type representing the type of data contained in a file.
    """
    CELLS = "cells", "The cells data type."
    CELL_SETS = "cell-sets", "The cell sets data type."
    EXPRESSION_MATRIX = "expression-matrix", "The gene expression matrix data type."
    MOLECULES = "molecules", "The molecules data type."
    NEIGHBORHOODS = "neighborhoods", "The spatial cell neighborhoods data type."
    RASTER = "raster", "The raster (i.e. imaging) data type."

class FileType(DocEnum):
    """
    An enum type representing the file format or schema to which a file conforms.
    """
    EXPRESSION_MATRIX_ZARR = "expression-matrix.zarr", "The Zarr-based expression matrix file type."
    CLUSTERS_JSON = "clusters.json", "The JSON-based expression matrix file type. Deprecated."
    GENES_JSON = "genes.json", "The JSON-based expression matrix file type. Deprecated."
    CELLS_JSON = "cells.json", "The JSON-based cells file type."
    MOLECULES_JSON = "molecules.json", "The JSON-based molecules file type."
    NEIGHBORHOODS_JSON = "neighborhoods.json", "The JSON-based neighborhoods file type."
    RASTER_JSON = "raster.json", "The JSON-based raster manifest file type."
    CELL_SETS_JSON = "cell-sets.json", "The JSON-based cell sets file type."
