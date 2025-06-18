#!/usr/bin/env python
# coding: utf-8

# # Vitessce Widget Tutorial

# ## Import dependencies
# 

# https://github.com/vitessce/vitessce/blob/main/examples/configs/src/view-configs/spatial-beta/kpmp.js

# In[2]:


from os.path import join


# In[1]:


import spatialdata
from spatialdata import SpatialData, to_polygons
from spatialdata.models import Image2DModel, Labels2DModel, TableModel
from spatialdata.transformations import (
    Affine,
    MapAxis,
    Scale,
    Sequence,
    Translation,
    get_transformation,
    set_transformation,
)
from anndata import read_zarr, AnnData
import tifffile
from tifffile import imread, TiffFile
import xml.etree.ElementTree
import io
import zarr



# In[2]:




# In[3]:


base_dir = join("..", "..", "..", "kpmp-f2f-march-2023", "S-1905-017737")


# In[ ]:


img_path = join(base_dir, "S-1905-017737_PAS_2of2_bf.ome.tif")
seg_path = join(base_dir, "S-1905-017737_PAS_2of2.ome.tif")

adata_paths = {
    "cortical_interstitia": join(base_dir, "Cortical Interstitium.adata.zarr"),
    "non_globally_sclerotic_glomeruli": join(base_dir, "Glomeruli.adata.zarr"),
    "globally_sclerotic_glomeruli": join(base_dir, "Globally Sclerotic Glomeruli.adata.zarr"),
    "tubules": join(base_dir, "Tubules with Area non infinity.adata.zarr"),
    "arteries_arterioles": None,
    "interstitialfibrosis_and_tubular_atrophy": join(base_dir, "IFTA.adata.zarr"),
    "peritubular_capillaries": join(base_dir, "Peritubular Capillaries renamed.adata.zarr"),
}


# In[ ]:


# The shape of the data should be c(z)yx for 2D (3D) images


# In[ ]:


img_arr = imread(img_path)


# In[ ]:


seg_store = imread(seg_path, aszarr=True)
seg_z = zarr.open(seg_store)


# In[ ]:


def clean_adata(adata):
    colnames = adata.obs.columns.tolist()
    adata.obs = adata.obs.rename(columns=dict([ (c, c.replace(" ", "_")) for c in colnames ]))
    # Many of the anndata objects contain un-helpful obs indices with entirely zero values like [0, 0, 0, 0]
    # TODO: determine whether using 1...N+1 in these cases is correct.
    if adata.obs.index.unique().shape[0] == 1:
        adata.obs.index = range(1, adata.obs.shape[0]+1)
    return adata


# In[ ]:


sdata = SpatialData(
    images={
        "image": Image2DModel.parse(img_arr, scale_factors=[2, 2, 2, 2, 2]),
    },
    tables={
        f"table_{obs_type}": TableModel.parse(clean_adata(read_zarr(adata_path)))
        for obs_type, adata_path in adata_paths.items()
        if adata_path is not None
    }
)

for i, obs_type in enumerate(adata_paths.keys()):
    sdata[f"labels_{obs_type}"] = Labels2DModel.parse(seg_z['1'][i, :, :], scale_factors=[2, 2, 2, 2])
    # Scale by a factor of 2 since we are using the second resolution from the original file.
    scale = Scale([2.0, 2.0], axes=("y","x"))
    set_transformation(sdata[f"labels_{obs_type}"], scale, to_coordinate_system="global")


# In[ ]:

for labels_key in sdata.labels.keys():
    shapes_key = "shapes" + labels_key[labels_key.index("_"):]
    sdata.shapes[shapes_key] = to_polygons(sdata.labels[labels_key])


sdata.write(join(base_dir, "sdata.zarr"), overwrite=True)
print(join(base_dir, "sdata.zarr"))
print("Done")

# In[ ]:





# In[ ]:




