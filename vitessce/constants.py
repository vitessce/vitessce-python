from enum import Enum

# Reference: https://stackoverflow.com/a/50473952


class DocEnum(Enum):
    def __new__(cls, value, doc):
        self = object.__new__(cls)
        self._value_ = value
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
    EMBEDDING_CELL_OPACITY = 'embeddingCellOpacity', "Opacity of cells in embedding (points in scatterplot view)"
    EMBEDDING_CELL_RADIUS = 'embeddingCellRadius', "Radius of cells in embedding (points in scatterplot view)"
    EMBEDDING_CELL_RADIUS_MODE = 'embeddingCellRadiusMode', "Radius mode of cells in embedding (points in scatterplot view) - auto or manual"
    EMBEDDING_CELL_OPACITY_MODE = 'embeddingCellOpacityMode', "Opacity mode of cells in embedding (points in scatterplot view) - auto or manual"
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
    CELL_FILTER = "cellFilter", "A subset of cells to include after filtering."
    CELL_HIGHLIGHT = "cellHighlight", "A subset of cells to highlight."
    CELL_SET_SELECTION = "cellSetSelection", "A subset of cell sets to select."
    CELL_SET_HIGHLIGHT = "cellSetHighlight", "A subset of cell sets to highlight."
    CELL_SET_COLOR = "cellSetColor", "A mapping from cell sets to colors."
    GENE_FILTER = "geneFilter", "A subset of genes to include after filtering."
    GENE_HIGHLIGHT = "geneHighlight", "A subset of genes to highlight."
    GENE_SELECTION = "geneSelection", "A subset of genes to select."
    GENE_EXPRESSION_COLORMAP = "geneExpressionColormap", "The colormap to use for the gene expression scale."
    GENE_EXPRESSION_COLORMAP_RANGE = "geneExpressionColormapRange", "The range of gene expression values to map."
    CELL_COLOR_ENCODING = "cellColorEncoding", "The color encoding to use for cell entities."
    SPATIAL_RASTER_LAYERS = 'spatialRasterLayers', "Layer definitions for the raster imagery in the spatial view."
    SPATIAL_CELLS_LAYER = 'spatialCellsLayer', "Layer definitions for the cells in the spatial view."
    SPATIAL_MOLECULES_LAYER = 'spatialMoleculesLayer', "Layer definitions for the molecules in the spatial view."
    SPATIAL_NEIGHBORHOODS_LAYER = 'spatialNeighborhoodsLayer', "Layer definitions for the neighborhoods in the spatial view."
    GENOMIC_ZOOM_X = "genomicZoomX", "The zoom level of a higlass view, X dimension."
    GENOMIC_ZOOM_Y = "genomicZoomY", "The zoom level of a higlass view, Y dimension."
    GENOMIC_TARGET_X = "genomicTargetX", "The x-coordinate of the center of a higlass view."
    GENOMIC_TARGET_Y = "genomicTargetY", "The y-coordinate of the center of a higlass view."
    ADDITIONAL_CELL_SETS = "additionalCellSets", "User-defined cell sets."


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
    GENOMIC_PROFILES = "genomicProfiles", "The higlass component can be used to visualize genome-wide ATAC-seq profiles."
    CELL_SET_SIZES = "cellSetSizes", "The cell set sizes component can display the quantities of cells in selected cell sets."
    GENES = "genes", "The gene list selector."
    CELL_SET_EXPRESSION = "cellSetExpression", "Expression levels are displayed by cell set"


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
    GENOMIC_PROFILES = "genomic-profiles", "The genomic profiles data type, used by HiGlass 1D quantitative tracks."


class FileType(DocEnum):
    """
    An enum type representing the file format or schema to which a file conforms.
    """
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
