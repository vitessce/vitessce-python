from os.path import join
from spatialdata import SpatialData, to_polygons
from spatialdata.models import Image2DModel, Labels2DModel, TableModel
from spatialdata.transformations import (
    Scale,
    set_transformation,
)
from anndata import read_zarr
from tifffile import imread
import zarr

base_dir = join("S-1905-017737")

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

# The shape of the data should be c(z)yx for 2D (3D) images
img_arr = imread(img_path)
seg_store = imread(seg_path, aszarr=True)
seg_z = zarr.open(seg_store)

def clean_adata(adata):
    colnames = adata.obs.columns.tolist()
    adata.obs = adata.obs.rename(columns=dict([ (c, c.replace(" ", "_")) for c in colnames ]))
    # Many of the anndata objects contain un-helpful obs indices with entirely zero values like [0, 0, 0, 0]
    # TODO: determine whether using 1...N+1 in these cases is correct.
    if adata.obs.index.unique().shape[0] == 1:
        adata.obs.index = range(1, adata.obs.shape[0]+1)
    return adata

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

for labels_key in sdata.labels.keys():
    shapes_key = "shapes" + labels_key[labels_key.index("_"):]
    sdata.shapes[shapes_key] = to_polygons(sdata.labels[labels_key])

# Convert boolean columns to set names
ptc_obs = sdata.tables['table_peritubular_capillaries'].obs
ptc_obs['PTC_in_Cortex'] = sdata.tables['table_peritubular_capillaries'].X[:,2]
ptc_obs['PTC_in_IFTA'] = sdata.tables['table_peritubular_capillaries'].X[:,3]

def row_to_cortex_name(row):
    if row['PTC_in_Cortex'] == 1.0:
        return "Inside Cortex"
    if row['PTC_in_Cortex'] == 0.0:
        return "Outside Cortex"

def row_to_ifta_name(row):
    if row['PTC_in_IFTA'] == 1.0:
        return "IFTA"
    if row['PTC_in_IFTA'] == 0.0:
        return "non-IFTA"

def row_to_cortical_ifta_name(row):
    if row['PTC_in_Cortex'] == 1.0 and row['PTC_in_IFTA'] == 1.0:
        return "Cortical IFTA"
    if row['PTC_in_Cortex'] == 1.0 and row['PTC_in_IFTA'] == 0.0:
        return "Cortical non-IFTA"
    if row['PTC_in_Cortex'] == 0.0:
        return "Outside Cortex"

ptc_obs['cortex_set'] = ptc_obs.apply(row_to_cortex_name, axis='columns')
ptc_obs['ifta_set'] = ptc_obs.apply(row_to_ifta_name, axis='columns')
ptc_obs['cortex_ifta_set'] = ptc_obs.apply(row_to_cortical_ifta_name, axis='columns')

sdata.tables['table_peritubular_capillaries'].obs = ptc_obs

sdata.write(join(base_dir, "sdata.zarr"), overwrite=True)
print(join(base_dir, "S-1905-017737.sdata.zarr"))
print("Done")
