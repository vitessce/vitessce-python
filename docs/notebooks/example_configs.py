spraggins = {
	"version": "0.1.0",
	"description": "High Bit Depth (uint16) Multiplex Immunofluorescence Imaging",
	"layers": [
		{
			"name": "raster",
			"type": "RASTER",
			"fileType": "raster.json",
			"url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/spraggins/spraggins.raster.json"
		}
	],
	"name": "Spraggins",
	"public": True,
	"staticLayout": [
		{
			"component": "spatial",
			"props": {
			"view": {
				"zoom": -6.5,
				"target": [
				20000,
				20000,
				0
				]
			}
			},
			"x": 0,
			"y": 0,
			"w": 9,
			"h": 2
		},
		{
			"component": "layerController",
			"x": 9,
			"y": 0,
			"w": 3,
			"h": 2
		}
	]
}
dries = {
    "name": "Dries",
	"version": "1.0.0",
	"description": "Giotto, a pipeline for integrative analysis and visualization of single-cell spatial transcriptomic data",
    "public": True,
    "datasets": [
        {
            "uid": 'dries-2019',
            "name": 'Dries 2019',
            "files": [
                {
                  "type": "cells",
                  "fileType": "cells.json",
                  "url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/dries/dries.cells.json"
                },
                {
                  "type": "cell-sets",
                  "fileType": "cell-sets.json",
                  "url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/dries/dries.cell-sets.json"
            }
            ]
        }
  	],
    "initStrategy": "auto",
    "coordinationSpace": {
      "embeddingType": {
        "TSNE": 't-SNE',
        "UMAP": 'UMAP',
      },
      "embeddingZoom": {
        "TSNE": 3,
        "UMAP": 3,
      },
      "spatialZoom": {
        "A": -4.4,
      },
      "spatialTargetX": {
        "A": 3800,
      },
      "spatialTargetY": {
        "A": 900,
      },
    },
    "layout": [
		{
		  "component": "description",
		  "props": {
			"description": "Giotto, a pipeline for integrative analysis and visualization of single-cell spatial transcriptomic data"
		  },
		  "x": 9,
		  "y": 0,
		  "w": 3,
		  "h": 2
		},
		{
		  "component": "status",
		  "x": 9,
		  "y": 2,
		  "w": 3,
		  "h": 2
		},
		{
		  "component": "cellSets",
		  "x": 9,
		  "y": 4,
		  "w": 3,
		  "h": 4
		},
		{
		  "component": "cellSetSizes",
		  "x": 5,
		  "y": 4,
		  "w": 4,
		  "h": 4
		},
		{
		  "component": "scatterplot",
		  "coordinationScopes": {
              "embeddingType": 'TSNE',
              "embeddingZoom": 'TSNE'
          },
          "x": 0,
		  "y": 2,
		  "w": 5,
		  "h": 4
		},
		{
		  "component": "spatial",
		  "props": {
			"cellRadius": 50
		  },
          "coordinationScopes": {
              "spatialZoom": 'A',
              "spatialTargetX": 'A',
              "spatialTargetY": 'A'
          },
		  "x": 5,
		  "y": 0,
		  "w": 4,
		  "h": 4
		},
		{
		  "component": "scatterplot",
		  "coordinationScopes": {
              "embeddingType": 'UMAP',
              "embeddingZoom": 'UMAP'
          },
		  "x": 0,
		  "y": 0,
		  "w": 5,
		  "h": 4
		}
  	]
}