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
