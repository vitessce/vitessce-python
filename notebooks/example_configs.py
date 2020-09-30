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
	"version": "0.1.0",
	"description": "Giotto, a pipeline for integrative analysis and visualization of single-cell spatial transcriptomic data",
    "layers": [
		{
		  "name": "cells",
		  "type": "CELLS",
		  "fileType": "cells.json",
		  "url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/dries/dries.cells.json"
		},
		{
		  "name": "cell-sets",
		  "type": "CELL-SETS",
		  "fileType": "cell-sets.json",
		  "url": "https://s3.amazonaws.com/vitessce-data/0.0.31/master_release/dries/dries.cell-sets.json"
		}
  	],
    "name": "Dries",
    "public": True,
    "staticLayout": [
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
		  "props": {
			"mapping": "t-SNE",
			"view": {
			  "zoom": 3,
			  "target": [
				0,
				0,
				0
			  ]
			}
		  },
		  "x": 0,
		  "y": 2,
		  "w": 5,
		  "h": 4
		},
		{
		  "component": "spatial",
		  "props": {
			"cellRadius": 50,
			"view": {
			  "zoom": -4.4,
			  "target": [
				3800,
				-900,
				0
			  ]
			}
		  },
		  "x": 5,
		  "y": 0,
		  "w": 4,
		  "h": 4
		},
		{
		  "component": "scatterplot",
		  "props": {
			"mapping": "UMAP",
			"view": {
			  "zoom": 3,
			  "target": [
				0,
				0,
				0
			  ]
			}
		  },
		  "x": 0,
		  "y": 0,
		  "w": 5,
		  "h": 4
		}
  	]
}