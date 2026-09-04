def get_next_scope_numeric(prev_scopes):
    next_scope_int = 0
    next_scope_str = None

    while True:
        next_scope_str = str(next_scope_int)
        if next_scope_str not in prev_scopes:
            break
        next_scope_int += 1
    return next_scope_str


def create_prefixed_get_next_scope_numeric(prefix):

    def inner_get_next_scope(prev_scopes):
        next_scope_int = 0
        next_scope_str = None

        while True:
            next_scope_str = f"{prefix}{next_scope_int}"
            if next_scope_str not in prev_scopes:
                break
            next_scope_int += 1
        return next_scope_str

    return inner_get_next_scope


def get_initial_coordination_scope_prefix(dataset_uid, data_type):
    return f"init_{dataset_uid}_{data_type}_"


def get_initial_coordination_scope_name(dataset_uid, data_type, i=None):
    prefix = get_initial_coordination_scope_prefix(dataset_uid, data_type)
    return f"{prefix}{0 if i is None else i}"


def make_ids_csv_data_url(ids, for_web_app=False):
    """
    Build a `data:` URL containing a small inline CSV with a single `id`
    column, given a list of observation IDs (e.g. segment IDs).

    Useful for defining a small, explicit `obsFeatureMatrix.csv` /
    `obsSets.csv`-style observation list without needing to host a
    separate CSV file.

    :param list ids: A list of observation IDs.
    :rtype: str
    :returns: A `data:text/csv,...` URL.
    """
    import csv
    import io
    from urllib.parse import quote

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["id"])
    writer.writerows([[i] for i in ids])
    encoded = quote(buf.getvalue())
    if for_web_app:
        encoded = quote(encoded, safe="")
    return f"data:text/csv,{encoded}"


def make_colors_csv_data_url(id_to_color, for_web_app=False):
    """
    Build a `data:` URL containing a small inline CSV with `id` and
    `color` columns, given a dict mapping observation ID to a color
    string (e.g. a hex color).

    Useful for defining explicit per-observation colors (`obsColors.csv`)
    without needing to host a separate CSV file.

    :param dict id_to_color: A dict mapping observation ID to color.
    :rtype: str
    :returns: A `data:text/csv,...` URL.
    """
    import csv
    import io
    from urllib.parse import quote

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["id", "color"])
    writer.writerows(id_to_color.items())
    encoded = quote(buf.getvalue())
    if for_web_app:
        encoded = quote(encoded, safe="")
    return f"data:text/csv,{encoded}"
