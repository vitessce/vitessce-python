from typing import Tuple, List, Optional

import os
from os.path import join
from bisect import bisect_left, bisect_right
import pandas as pd
import numpy as np

from spatialdata import get_element_annotators
from spatialdata.models import PointsModel
import dask.dataframe as dd


MORTON_CODE_NUM_BITS = 32  # Resulting morton codes will be stored as uint32.
MORTON_CODE_VALUE_MIN = 0
MORTON_CODE_VALUE_MAX = 2**(MORTON_CODE_NUM_BITS / 2) - 1

# In the first few rows (up to four), we set the morton code
# value to a sentinel value for the rows that contain extrema,
# x_min, x_max, y_min, and/or y_max data extent.
# It is possible for the same row to contain extrema along more than one axis,
# in which case we will only use fewer than four rows for the extrema.
MORTON_CODE_EXTREME_VALUE_INDICATOR = 0

# --------------------------
# Functions for computing Morton codes for SpatialData points (2D).
# --------------------------


def norm_series_to_uint(series, v_min, v_max):
    """
    Scale numeric Series (int or float) to integer grid [0, 2^bits-1], handling NaNs.
    """
    # Cast to float64
    series_f64 = series.astype("float64")
    # Normalize the array values to be between 0.0 and 1.0
    norm_series_f64 = (series_f64 - v_min) / (v_max - v_min)
    # Clip to ensure no values are outside 0/1 range
    clipped_norm_series_f64 = np.clip(norm_series_f64, 0.0, 1.0)
    # Multiply by the morton code max-value to scale from [0,1] to [0,65535]
    out = (clipped_norm_series_f64 * MORTON_CODE_VALUE_MAX).astype(np.uint32)
    # Set NaNs to 0.
    out = out.fillna(0)
    return out


def norm_ddf_to_uint(ddf):
    [x_min, x_max, y_min, y_max] = [ddf["x"].min().compute(), ddf["x"].max().compute(), ddf["y"].min().compute(), ddf["y"].max().compute()]
    ddf["x_uint"] = norm_series_to_uint(ddf["x"], x_min, x_max)
    ddf["y_uint"] = norm_series_to_uint(ddf["y"], y_min, y_max)
    return ddf


def _part1by1_16(x):
    """
    Spread each 16-bit value into 32 bits by inserting zeros between bits.
    Input:  uint32 array (values must fit in 16 bits)
    Output: uint32 array (bit-spread)
    """

    assert x.dtype.name == 'uint32'

    # Mask away any bits above 16 (just in case input wasn't clean).
    x = x & np.uint32(0x0000FFFF)

    # First spread: shift left by 8 bits, OR with original, then mask.
    # After this, groups of 8 bits are separated by 8 zeros.
    x = (x | np.left_shift(x, 8)) & np.uint32(0x00FF00FF)

    # Spread further: now groups of 4 bits separated by 4 zeros.
    x = (x | np.left_shift(x, 4)) & np.uint32(0x0F0F0F0F)

    # Spread further: groups of 2 bits separated by 2 zeros.
    x = (x | np.left_shift(x, 2)) & np.uint32(0x33333333)

    # Final spread: single bits separated by a zero bit.
    # Now each original bit is in every other position (positions 0,2,4,...).
    x = (x | np.left_shift(x, 1)) & np.uint32(0x55555555)

    return x


def _part1by1_32(x):
    """
    Spread each 32-bit value into 64 bits by inserting zeros between bits.
    Input:  uint64 array (values must fit in 32 bits)
    Output: uint64 array (bit-spread)
    """

    assert x.dtype.name == 'uint64'

    # Mask away any bits above 32 (safety).
    x = x.astype(np.uint64) & np.uint64(0x00000000FFFFFFFF)

    # First spread: separate into 16-bit chunks spaced out.
    x = (x | np.left_shift(x, 16)) & np.uint64(0x0000FFFF0000FFFF)

    # Spread further: each 8-bit chunk separated.
    x = (x | np.left_shift(x, 8)) & np.uint64(0x00FF00FF00FF00FF)

    # Spread further: each 4-bit nibble separated.
    x = (x | np.left_shift(x, 4)) & np.uint64(0x0F0F0F0F0F0F0F0F)

    # Spread further: 2-bit groups separated.
    x = (x | np.left_shift(x, 2)) & np.uint64(0x3333333333333333)

    # Final spread: single bits separated by zeros.
    # Now each original bit occupies every other position (0,2,4,...).
    x = (x | np.left_shift(x, 1)) & np.uint64(0x5555555555555555)

    return x


def morton_interleave(ddf):
    """
    Vectorized Morton interleave for integer arrays xi, yi
    already scaled to [0, 2^bits - 1].
    Returns Morton codes as uint32 (if bits<=16) or uint64 (if bits<=32).
    """

    xi = ddf["x_uint"]
    yi = ddf["y_uint"]

    # Spread x and y bits into even (x) and odd (y) positions.
    xs = _part1by1_16(xi)
    ys = _part1by1_16(yi)

    # Interleave: shift y bits left by 1 so they go into odd positions,
    # then OR with x bits in even positions.
    code = np.left_shift(ys.astype(np.uint64), 1) | xs.astype(np.uint64)

    # Fits in 32 bits since we only had 16+16 input bits.
    return code.astype(np.uint32)


def sdata_morton_sort_points(sdata, element):
    ddf = sdata.points[element]

    # Add normalized integer columns and morton codes using the full dataset bounding box.
    ddf = norm_ddf_to_uint(ddf)
    ddf["morton_code_2d"] = morton_interleave(ddf)

    # Identify the extreme coordinate values.
    x_min_val = ddf["x"].min().compute()
    x_max_val = ddf["x"].max().compute()
    y_min_val = ddf["y"].min().compute()
    y_max_val = ddf["y"].max().compute()

    # Build a boolean mask for the extreme rows using coordinate value equality.
    is_extreme = (
        (ddf["x"] == x_min_val) | (ddf["x"] == x_max_val)
        | (ddf["y"] == y_min_val) | (ddf["y"] == y_max_val)
    )

    extreme_pdf = ddf.loc[is_extreme].compute().reset_index(drop=True)

    # Checking based on the row indices would be cleaner, but dask does not
    # guarantee unique indices across partitions.
    # Reference: https://docs.dask.org/en/stable/generated/dask.dataframe.DataFrame.reset_index.html
    # Instead, our `is_extreme` uses the min/max values themselves to identify the matching rows,
    # which can potentially return more than four rows if multiple rows match these extreme values.
    # TODO: in that case, reorder the extreme rows so the first four form the bounding box,
    # and any additional rows should be moved into `rest_df` before the morton-based sorting occurs.
    assert extreme_pdf.shape[0] <= 4

    # Mark sentinel rows with morton_code_2d = MORTON_CODE_EXTREME_VALUE_INDICATOR so the reader can identify them without
    # relying on a hard-coded count or fixed row positions.
    extreme_pdf["morton_code_2d"] = np.uint32(MORTON_CODE_EXTREME_VALUE_INDICATOR)

    # Exclude the extreme rows from the rest.
    rest_ddf = ddf.loc[~is_extreme]

    if "z" in ddf.columns:
        num_unique_z = ddf["z"].unique().shape[0].compute()
        if num_unique_z < 100:
            # Heuristic for interpreting the 3D data as 2.5D
            # Reference: https://github.com/scverse/spatialdata/issues/961
            sorted_rest = rest_ddf.sort_values(by=["z", "morton_code_2d"], ascending=True)
        else:
            # TODO: include z as a dimension in the morton code in the 3D case?
            sorted_rest = rest_ddf.sort_values(by="morton_code_2d", ascending=True)
    else:
        sorted_rest = rest_ddf.sort_values(by="morton_code_2d", ascending=True)

    # Prepend the sentinel rows (morton_code_2d == MORTON_CODE_EXTREME_VALUE_INDICATOR) before the morton-sorted data.
    sentinel_ddf = dd.from_pandas(extreme_pdf, npartitions=1)
    result_ddf = dd.concat([sentinel_ddf, sorted_rest], ignore_index=True, interleave_partitions=True)

    sdata.points[element] = PointsModel.parse(result_ddf)

    # annotating_tables = get_element_annotators(sdata, element)

    # TODO: Sort any annotating table(s) as well.

    return sdata


def sdata_morton_query_rect_aux(sdata, element, orig_rect):
    # orig_rect = [[50, 50], [100, 150]] # [[x0, y0], [x1, y1]]
    # norm_rect = [
    #    orig_coord_to_norm_coord(orig_rect[0], orig_x_min=0, orig_x_max=100, orig_y_min=0, orig_y_max=200),
    #    orig_coord_to_norm_coord(orig_rect[1], orig_x_min=0, orig_x_max=100, orig_y_min=0, orig_y_max=200)
    # ]

    sorted_ddf = sdata.points[element]

    # Sentinel rows are identified by morton_code_2d == 0 and always precede the z-order data.
    # There are at most 4 (one per axis extreme), fewer when a single point covers multiple extremes.
    head = sorted_ddf.head(4)
    sentinel_rows = head[head["morton_code_2d"] == MORTON_CODE_EXTREME_VALUE_INDICATOR]
    x_min = float(sentinel_rows["x"].min())
    x_max = float(sentinel_rows["x"].max())
    y_min = float(sentinel_rows["y"].min())
    y_max = float(sentinel_rows["y"].max())

    norm_rect = [
        orig_coord_to_norm_coord(orig_rect[0], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max),
        orig_coord_to_norm_coord(orig_rect[1], orig_x_min=x_min, orig_x_max=x_max, orig_y_min=y_min, orig_y_max=y_max)
    ]

    # Get a list of morton code intervals that cover this rectangle region
    # [ (morton_start, morton_end), ... ]
    morton_intervals = zcover_rectangle(
        rx0=norm_rect[0][0], ry0=norm_rect[0][1],
        rx1=norm_rect[1][0], ry1=norm_rect[1][1],
        bits=16,
        stop_level=None,
        merge=True,
    )

    return morton_intervals


def sdata_morton_query_rect(sdata, element, orig_rect):
    sorted_ddf = sdata.points[element]

    # TODO: generalize to 3D morton codes

    morton_intervals = sdata_morton_query_rect_aux(sdata, element, orig_rect)

    # Get morton code column as a list of integers.
    morton_sorted = sorted_ddf["morton_code_2d"].compute().values.tolist()

    # Count leading sentinel rows (morton_code_2d == 0); stop at the first non-zero code.
    sentinel_count = 0
    for code in morton_sorted[:4]:
        if code == MORTON_CODE_EXTREME_VALUE_INDICATOR:
            sentinel_count += 1
        else:
            break

    # Get a list of row ranges that match the morton intervals.
    # (This uses binary searches internally to find the matching row indices).
    # [ (row_start, row_end), ... ]
    matching_row_ranges = zquery_rows(morton_sorted[sentinel_count:], morton_intervals, merge=True)

    # Offset ranges to account for the sentinel rows at the start of the table.
    matching_row_ranges = [(i + sentinel_count, j + sentinel_count) for i, j in matching_row_ranges]

    return matching_row_ranges


def sdata_morton_query_rect_debug(sdata, element, orig_rect):
    # This is the same as the above sdata_morton_query_rect function,
    # but it also returns the list of row indices that were checked
    # during the binary searches.
    sorted_ddf = sdata.points[element]
    morton_intervals = sdata_morton_query_rect_aux(sdata, element, orig_rect)
    morton_sorted = sorted_ddf["morton_code_2d"].compute().values.tolist()
    sentinel_count = 0
    for code in morton_sorted[:4]:
        if code == MORTON_CODE_EXTREME_VALUE_INDICATOR:
            sentinel_count += 1
        else:
            break
    matching_row_ranges, rows_checked = zquery_rows_aux(morton_sorted[sentinel_count:], morton_intervals, merge=True)
    matching_row_ranges = [(i + sentinel_count, j + sentinel_count) for i, j in matching_row_ranges]
    return matching_row_ranges, rows_checked

# --------------------------
# Functions for rectangle queries.
# --------------------------

# Convert a coordinate from the normalized [0, 65535] space to the original space.


def norm_coord_to_orig_coord(norm_coord, orig_x_min, orig_x_max, orig_y_min, orig_y_max):
    [norm_x, norm_y] = norm_coord
    orig_x_range = orig_x_max - orig_x_min
    orig_y_range = orig_y_max - orig_y_min
    return [
        (orig_x_min + (norm_x / MORTON_CODE_VALUE_MAX) * orig_x_range),
        (orig_y_min + (norm_y / MORTON_CODE_VALUE_MAX) * orig_y_range),
    ]

# Convert a coordinate from the original space to the [0, 65535] normalized space.


def orig_coord_to_norm_coord(orig_coord, orig_x_min, orig_x_max, orig_y_min, orig_y_max):
    [orig_x, orig_y] = orig_coord
    orig_x_range = orig_x_max - orig_x_min
    orig_y_range = orig_y_max - orig_y_min
    return [
        np.float64(((orig_x - orig_x_min) / orig_x_range) * MORTON_CODE_VALUE_MAX).astype(np.uint32),
        np.float64(((orig_y - orig_y_min) / orig_y_range) * MORTON_CODE_VALUE_MAX).astype(np.uint32),
    ]

# --------------------------
# Quadtree / Z-interval helpers
# --------------------------


def intersects(ax0: int, ay0: int, ax1: int, ay1: int,
               bx0: int, by0: int, bx1: int, by1: int) -> bool:
    """Axis-aligned box intersection (inclusive integer bounds)."""
    return not (ax1 < bx0 or bx1 < ax0 or ay1 < by0 or by1 < ay0)


def contained(ix0: int, iy0: int, ix1: int, iy1: int,
              ox0: int, oy0: int, ox1: int, oy1: int) -> bool:
    """Is inner box entirely inside outer box? (inclusive integer bounds)"""
    return (ox0 <= ix0 <= ix1 <= ox1) and (oy0 <= iy0 <= iy1 <= oy1)


def point_inside(x: int, y: int, rx0: int, ry0: int, rx1: int, ry1: int) -> bool:
    return (rx0 <= x <= rx1) and (ry0 <= y <= ry1)


def cell_range(prefix: int, level: int, bits: int) -> Tuple[int, int]:
    """
    All Morton codes in a quadtree cell share the same prefix (2*level bits).
    Fill the remaining lower bits with 0s (lo) or 1s (hi).
    """
    shift = 2 * (bits - level)
    lo = prefix << shift
    hi = ((prefix + 1) << shift) - 1
    return lo, hi


def merge_adjacent(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping or directly adjacent intervals."""
    if not intervals:
        return []
    intervals.sort(key=lambda t: t[0])
    merged = [intervals[0]]
    for lo, hi in intervals[1:]:
        mlo, mhi = merged[-1]
        if lo <= mhi + 1:
            merged[-1] = (mlo, max(mhi, hi))
        else:
            merged.append((lo, hi))
    return merged

# --------------------------
# Rectangle -> list of Morton intervals
# --------------------------


def zcover_rectangle(rx0: int, ry0: int, rx1: int, ry1: int, bits: int, stop_level: Optional[int] = None, merge: bool = True) -> List[Tuple[int, int]]:
    """
    Compute a (near-)minimal set of Morton code ranges covering the rectangle
    [rx0..rx1] x [ry0..ry1] on an integer grid [0..2^bits-1]^2.

    - If stop_level is None: exact cover (descend to exact containment).
    - If stop_level is set (0..bits): stop descending at that level, adding
      partially-overlapping cells as whole ranges (superset cover).
    """
    if not (0 <= rx0 <= rx1 <= (1 << bits) - 1 and 0 <= ry0 <= ry1 <= (1 << bits) - 1):
        raise ValueError("Rectangle out of bounds for given bits.")

    intervals: List[Tuple[int, int]] = []

    # stack entries: (prefix, level, xmin, ymin, xmax, ymax)
    stack = [(0, 0, 0, 0, (1 << bits) - 1, (1 << bits) - 1)]

    while stack:
        prefix, level, xmin, ymin, xmax, ymax = stack.pop()

        if not intersects(xmin, ymin, xmax, ymax, rx0, ry0, rx1, ry1):
            continue

        # If we stop at this level for a loose cover, add full cell range.
        if stop_level is not None and level == stop_level:
            intervals.append(cell_range(prefix, level, bits))
            continue

        # Fully contained: add full cell range.
        if contained(xmin, ymin, xmax, ymax, rx0, ry0, rx1, ry1):
            intervals.append(cell_range(prefix, level, bits))
            continue

        # Leaf cell: single lattice point (only happens when level==bits)
        if level == bits:
            if point_inside(xmin, ymin, rx0, ry0, rx1, ry1):
                intervals.append(cell_range(prefix, level, bits))
            continue

        # Otherwise, split into 4 children (Morton order: 00,01,10,11)
        midx = (xmin + xmax) // 2
        midy = (ymin + ymax) // 2

        # q0: (x<=midx, y<=midy) -> child code 0b00
        stack.append(((prefix << 2) | 0,
                      level + 1,
                      xmin, ymin, midx, midy))
        # q1: (x>midx, y<=midy)  -> child code 0b01
        stack.append(((prefix << 2) | 1,
                      level + 1,
                      midx + 1, ymin, xmax, midy))
        # q2: (x<=midx, y>midy)  -> child code 0b10
        stack.append(((prefix << 2) | 2,
                      level + 1,
                      xmin, midy + 1, midx, ymax))
        # q3: (x>midx, y>midy)   -> child code 0b11
        stack.append(((prefix << 2) | 3,
                      level + 1,
                      midx + 1, midy + 1, xmax, ymax))

    return merge_adjacent(intervals) if merge else intervals


# --------------------------
# Morton intervals -> row ranges in a Morton-sorted column
# --------------------------

def zquery_rows_aux(morton_sorted: List[int], intervals: List[Tuple[int, int]], merge: bool = True) -> Tuple[List[Tuple[int, int]], List[int]]:
    """
    For each Z-interval [zlo, zhi], binary-search in the sorted Morton column
    and return row index half-open ranges [i, j) to scan.
    """

    # Keep track of which keys were looked at during the binary searches.
    # This is used for analysis / debugging, for instance, to enable
    # evaluating how many HTTP requests would be needed in network-based case
    # (which will also depend on Arrow row group size).
    recorded_keys = []

    def record_key_check(k: int) -> int:
        # TODO: Does recorded_keys need to be marked as a global here?
        recorded_keys.append(k)
        return k

    ranges: List[Tuple[int, int]] = []
    # TODO: can these multiple binary searches be optimized?
    # Since we are doing many searches in the same array, and in each search we learn where more elements are located.
    for zlo, zhi in intervals:
        i = bisect_left(morton_sorted, zlo, key=record_key_check)
        # TODO: use lo=i in bisect_right to limit the search range?
        # TODO: can the second binary search be further optimized since we just did a binary search via bisect_left?
        j = bisect_right(morton_sorted, zhi, key=record_key_check)
        if i < j:
            ranges.append((i, j))

    result = merge_adjacent(ranges) if merge else ranges
    return result, recorded_keys


def zquery_rows(morton_sorted: List[int], intervals: List[Tuple[int, int]], merge: bool = True) -> List[Tuple[int, int]]:
    """
    For each Z-interval [zlo, zhi], binary-search in the sorted Morton column
    and return row index half-open ranges [i, j) to scan.
    """
    return zquery_rows_aux(morton_sorted, intervals, merge=merge)[0]


def row_ranges_to_row_indices(intervals: List[Tuple[int, int]]) -> List[int]:
    """
    Convert row ranges [i, j) to a list of row indices.
    Then, can index into pandas DataFrame using df.iloc[indices, :]
    """
    indices: List[int] = []
    for i, j in intervals:
        indices.extend(list(range(i, j)))
    return indices


# More helper functions.
def sdata_points_process_columns(sdata, element, var_name_col=None, table_name=None) -> dd.DataFrame:
    ddf = sdata.points[element]

    if var_name_col is None:
        # We can try to get it from the spatialdata_attrs metadata.
        var_name_col = sdata.points[element].attrs["spatialdata_attrs"].get("feature_key")

    # Appending codes for dictionary-encoded feature_name column.
    if table_name is None and var_name_col is not None:
        annotating_tables = get_element_annotators(sdata, element)
        if len(annotating_tables) == 1:
            table_name = annotating_tables[0]
        elif len(annotating_tables) == 0:
            raise ValueError(f"No annotating table found for Points element {element}, please specify table_name explicitly.")
        else:
            raise ValueError(f"Multiple annotating tables found for Points element {element}, please specify table_name explicitly.")

    if var_name_col is not None:
        var_df = sdata.tables[table_name].var
        var_index = var_df.index.values.tolist()

        def try_index(gene_name):
            try:
                return var_index.index(gene_name)
            except BaseException:
                return -1
        ddf[f"{var_name_col}_codes"] = ddf[var_name_col].apply(try_index).astype('int32')

    # Identify dictionary-encoded columns (categorical/string)
    orig_columns = ddf.columns.tolist()
    dict_encoded_cols = [col for col in orig_columns if pd.api.types.is_categorical_dtype(ddf[col].dtype) or pd.api.types.is_string_dtype(ddf[col].dtype)]

    # Dictionary-encoded columns (i.e., categorical and string) must be stored as the rightmost columns of the dataframe.
    ordered_columns = sorted(orig_columns, key=lambda colname: orig_columns.index(colname) if colname not in dict_encoded_cols else len(orig_columns))

    # Reorder the columns of the dataframe
    ddf = ddf[ordered_columns]

    return ddf


def sdata_points_modify_row_group_size(sdata, element, row_group_size: int = 50_000):
    import pyarrow.parquet as pq

    sdata_path = sdata.path
    # TODO: error if no path

    # List the parts of the parquet file.
    parquet_path = join(sdata_path, "points", element, "points.parquet")

    # Read the number of "part.*.parquet" files on disk.
    part_files = [f for f in os.listdir(parquet_path) if f.startswith("part.") and f.endswith(".parquet")]
    num_parts = len(part_files)

    # Update the row group size in each .parquet file part.
    for i in range(num_parts):
        part_path = join(parquet_path, f"part.{i}.parquet")
        table_read = pq.read_table(part_path)

        # Write the table to a new Parquet file with the desired row group size.
        pq.write_table(table_read, part_path, row_group_size=row_group_size)
