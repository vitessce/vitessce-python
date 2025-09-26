import pytest
from os.path import join
import numpy as np

from spatialdata import read_zarr

from vitessce.data_utils.spatialdata_points_zorder import (
    # Function for computing codes and sorting
    sdata_morton_sort_points,
    # Functions for querying
    sdata_morton_query_rect_debug,
    row_ranges_to_row_indices,
    orig_coord_to_norm_coord,
    MORTON_CODE_VALUE_MAX,
)

def is_sorted(l):
    return all(l[i] <= l[i + 1] for i in range(len(l) - 1))


@pytest.fixture
def sdata_with_points():
    data_dir = join("docs", "notebooks", "data")
    spatialdata_filepath = join(data_dir, "xenium_rep1_io.spatialdata.zarr")

    sdata = read_zarr(spatialdata_filepath)
    return sdata


def test_zorder_sorting(sdata_with_points):
    sdata = sdata_with_points

    sdata_morton_sort_points(sdata, "transcripts")

    # Check that the morton codes are sorted
    sorted_ddf = sdata.points["transcripts"]
    morton_sorted = sorted_ddf["morton_code_2d"].compute().values.tolist()

    assert is_sorted(morton_sorted)


def norm_value_to_uint(value, v_min, v_max):
    """
    Scale numeric value (int or float) to integer [0, 2^bits-1].
    """
    # Cast to float64
    value_f64 = value.astype("float64")
    # Normalize the array values to be between 0.0 and 1.0
    norm_value_f64 = (value_f64 - v_min) / (v_max - v_min)
    # Clip to ensure no values are outside 0/1 range
    clipped_norm_value_f64 = np.clip(norm_value_f64, 0.0, 1.0)
    # Multiply by the morton code max-value to scale from [0,1] to [0,65535]
    out = (clipped_norm_value_f64 * MORTON_CODE_VALUE_MAX).astype(np.uint32)
    return out


def test_zorder_query(sdata_with_points):
    sdata = sdata_with_points

    sdata_morton_sort_points(sdata, "transcripts")

    # Query a rectangle that should return some points
    orig_rect = [[50.0, 50.0], [100.0, 150.0]] # x0, y0, x1, y1
    matching_row_ranges, rows_checked = sdata_morton_query_rect_debug(sdata, "transcripts", orig_rect)
    rect_row_indices = row_ranges_to_row_indices(matching_row_ranges)

    # Cannot use df.iloc on a dask dataframe, so convert it to pandas first
    ddf = sdata.points["transcripts"]
    df = ddf.compute()
    df = df.reset_index(drop=True)
    estimated_row_indices = df.iloc[rect_row_indices].index.tolist()

    assert df.shape[0] == 42638083

    # Do the same query the "dumb" way, by checking all points

    # We need an epsilon for the "dumb" query since the normalization
    # introduces rounding issues. We can instead verify that a slightly
    # smaller rectangle is fully contained in the morton code query
    # estimated results.
    EXACT_BOUNDARY_EPSILON = 1

    in_rect = (
        (df["x"] >= orig_rect[0][0] + EXACT_BOUNDARY_EPSILON)
        & (df["x"] <= orig_rect[1][0] - EXACT_BOUNDARY_EPSILON)
        & (df["y"] >= orig_rect[0][1] + EXACT_BOUNDARY_EPSILON)
        & (df["y"] <= orig_rect[1][1] - EXACT_BOUNDARY_EPSILON)
    )
    dumb_df_subset = df.loc[in_rect]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices = dumb_df_subset.index.tolist()

    # Check that the estimated rows 100% contain the exact rows.
    # A.issubset(B) checks that all elements of A are in B ("A is a subset of B").
    assert set(exact_row_indices).issubset(set(estimated_row_indices))
    assert len(exact_row_indices) == 552
    assert len(estimated_row_indices) <= 631

    # Check that the number of rows checked is less than the total number of points
    assert len(rows_checked) <= 85374
    assert len(matching_row_ranges) == 24 # Kind of an implementation detail.

    # Do a second check, this time against x_uint/y_uint (the normalized coordinates)
    # TODO: does this ensure that estimated == exact?

    bounding_box = ddf.attrs["bounding_box"]
    x_min = bounding_box["x_min"]
    x_max = bounding_box["x_max"]
    y_min = bounding_box["y_min"]
    y_max = bounding_box["y_max"]
    norm_rect = [
        orig_coord_to_norm_coord(orig_rect[0], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max),
        orig_coord_to_norm_coord(orig_rect[1], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max)
    ]

    in_rect_norm = (
        (df["x_uint"] >= norm_rect[0][0])
        & (df["x_uint"] <= norm_rect[1][0])
        & (df["y_uint"] >= norm_rect[0][1])
        & (df["y_uint"] <= norm_rect[1][1])
    )
    dumb_df_subset_norm = df.loc[in_rect_norm]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices_norm = dumb_df_subset_norm.index.tolist()

    # A.issubset(B)
    # True if A is a subset of B and False otherwise.
    assert set(exact_row_indices_norm).issubset(set(estimated_row_indices))

    assert len(exact_row_indices_norm) == 618
    assert len(estimated_row_indices) <= 631

    
    # ========= Another query ==========
    orig_rect = [[500.0, 500.0], [600.0, 600.0]] # x0, y0, x1, y1

    # Query using z-order
    matching_row_ranges, rows_checked = sdata_morton_query_rect_debug(sdata, "transcripts", orig_rect)
    rect_row_indices = row_ranges_to_row_indices(matching_row_ranges)
    estimated_row_indices = df.iloc[rect_row_indices].index.tolist()

    # Do the same query the "dumb" way, by checking all points
    in_rect = (
        (df["x"] >= orig_rect[0][0] + EXACT_BOUNDARY_EPSILON)
        & (df["x"] <= orig_rect[1][0] - EXACT_BOUNDARY_EPSILON)
        & (df["y"] >= orig_rect[0][1] + EXACT_BOUNDARY_EPSILON)
        & (df["y"] <= orig_rect[1][1] - EXACT_BOUNDARY_EPSILON)
    )
    dumb_df_subset = df.loc[in_rect]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices = dumb_df_subset.index.tolist()

    # Check that the estimated rows 100% contain the exact rows.
    # A.issubset(B) checks that all elements of A are in B ("A is a subset of B").
    assert set(exact_row_indices).issubset(set(estimated_row_indices))
    assert len(exact_row_indices) == 16678
    assert len(estimated_row_indices) <= 17681

    # Check that the number of rows checked is less than the total number of points
    assert len(rows_checked) <= 124661
    assert len(matching_row_ranges) == 176 # Kind of an implementation detail.

    # Do the same query the "dumb" way, by checking all points
    in_rect = (
        (df["x"] >= orig_rect[0][0] + EXACT_BOUNDARY_EPSILON)
        & (df["x"] <= orig_rect[1][0] - EXACT_BOUNDARY_EPSILON)
        & (df["y"] >= orig_rect[0][1] + EXACT_BOUNDARY_EPSILON)
        & (df["y"] <= orig_rect[1][1] - EXACT_BOUNDARY_EPSILON)
    )
    dumb_df_subset = df.loc[in_rect]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices = dumb_df_subset.index.tolist()

    # Query 2: Do a second check, this time against x_uint/y_uint (the normalized coordinates)
    norm_rect = [
        orig_coord_to_norm_coord(orig_rect[0], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max),
        orig_coord_to_norm_coord(orig_rect[1], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max)
    ]

    in_rect_norm = (
        (df["x_uint"] >= norm_rect[0][0])
        & (df["x_uint"] <= norm_rect[1][0])
        & (df["y_uint"] >= norm_rect[0][1])
        & (df["y_uint"] <= norm_rect[1][1])
    )
    dumb_df_subset_norm = df.loc[in_rect_norm]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices_norm = dumb_df_subset_norm.index.tolist()

    # A.issubset(B)
    # True if A is a subset of B and False otherwise.
    assert set(exact_row_indices_norm).issubset(set(estimated_row_indices))

    # Check that the estimated rows contain all of the exact rows.
    assert len(exact_row_indices_norm) == 17590
    assert len(estimated_row_indices) <= 17681
    

