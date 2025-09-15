import pytest
from os.path import join

from spatialdata import read_zarr

from vitessce.data_utils.spatialdata_points_zorder import (
    # Function for computing codes and sorting
    sdata_morton_sort_points,
    # Functions for querying
    sdata_morton_query_rect_debug,
    row_ranges_to_row_indices,
    orig_coord_to_norm_coord,
)

def is_sorted(l):
    return all(l[i] <= l[i + 1] for i in range(len(l) - 1))

def get_sdata():
    data_dir = join("docs", "notebooks", "data")
    spatialdata_filepath = join(data_dir, "xenium_rep1_io.spatialdata.zarr")

    sdata = read_zarr(spatialdata_filepath)
    return sdata

@pytest.mark.skip(reason="Temporarily disable")
def test_zorder_sorting():
    # TODO: use fixture here
    sdata = get_sdata()

    sdata_morton_sort_points(sdata, "transcripts")

    # Check that the morton codes are sorted
    sorted_ddf = sdata.points["transcripts"]
    morton_sorted = sorted_ddf["morton_code_2d"].compute().values.tolist()

    assert is_sorted(morton_sorted)


def test_zorder_query():
    sdata = get_sdata()

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
    in_rect = (
        (df["x"] >= orig_rect[0][0])
        & (df["x"] <= orig_rect[1][0])
        & (df["y"] >= orig_rect[0][1])
        & (df["y"] <= orig_rect[1][1])
    )
    dumb_df_subset = df.loc[in_rect]
    # Get the row indices of the points in the rectangle
    # (these are the indices in the original dataframe)
    exact_row_indices = dumb_df_subset.index.tolist()

    # Check that the estimated rows 100% contain the exact rows.
    # A.issubset(B) checks that all elements of A are in B ("A is a subset of B").
    assert set(exact_row_indices).issubset(set(estimated_row_indices))
    assert len(exact_row_indices) == 614
    assert len(estimated_row_indices) <= 631

    # Check that the number of rows checked is less than the total number of points
    assert len(rows_checked) <= 45237
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
    assert set(exact_row_indices_norm).issubset(set(estimated_row_indices))
    assert len(exact_row_indices_norm) == 617
    assert len(estimated_row_indices) <= 631
    


    """
    # ========= Another query ==========
    orig_rect = [[500, 500], [600, 600]] # x0, y0, x1, y1

    # Query using z-order
    matching_row_ranges, rows_checked = sdata_morton_query_rect_debug(sdata, "transcripts", orig_rect)
    rect_row_indices = row_ranges_to_row_indices(matching_row_ranges)
    estimated_row_indices = df.iloc[rect_row_indices].index.tolist()

    # Query the "dumb" way
    in_rect = (
        (df["x"] >= orig_rect[0][0])
        & (df["x"] <= orig_rect[1][0])
        & (df["y"] >= orig_rect[0][1])
        & (df["y"] <= orig_rect[1][1])
    )
    dumb_df_subset = df.loc[in_rect]
    exact_row_indices = dumb_df_subset.index.tolist()

    diff_rows = set(estimated_row_indices) - set(exact_row_indices)
    # print("Rows in estimated but not exact:", diff_rows)
    print(df.iloc[list(diff_rows)])
    raise NotImplementedError("Debugging")

    # Check that the estimated rows contain all of the exact rows.
    assert len(set(exact_row_indices).intersection(set(estimated_row_indices))) == 0
    assert len(exact_row_indices) <= 1123 # TODO: update
    assert len(estimated_row_indices) <= 1163 # TODO: update
    
    """







    
