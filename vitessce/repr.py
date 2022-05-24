import inspect


def make_repr(init_locals, class_def=None):
    '''
    >>> from .wrappers import MultiImageWrapper
    >>> orig = MultiImageWrapper('IMAGE_WRAPPERS', foo='bar')
    >>> orig_repr = repr(orig)
    >>> print(orig_repr)
    MultiImageWrapper(image_wrappers='IMAGE_WRAPPERS', foo='bar')
    >>> evalled = eval(orig_repr)
    >>> assert orig_repr == repr(evalled)
    '''
    # Get the class definition from locals.
    clazz = None
    if '__class__' in init_locals:
        # Requires superclass to be initialized.
        clazz = init_locals.pop('__class__')
    elif 'self' in init_locals and hasattr(init_locals['self'], '__class__'):
        clazz = init_locals["self"].__class__  # pragma: no cover
    elif class_def is not None:
        clazz = class_def
    else:
        raise ValueError("make_repr could not locate the class definition")  # pragma: no cover

    # Remove self from locals.
    if 'self' in init_locals:
        del init_locals['self']

    # Get the class name.
    class_name = clazz.__name__

    # Remove redundant constructor parameters (when the value equals the default value).
    for k, v in inspect.signature(clazz).parameters.items():
        if k not in init_locals:
            continue

        try:
            is_default = init_locals[k] == v.default
        except NotImplementedError:
            # Equality comparison may not be implemented for the value object.
            is_default = False

        try:
            if is_default:
                del init_locals[k]
        except ValueError:
            # TODO: Is this also expected?
            # ValueError('The truth value of a DataFrame is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().')
            pass

    # Convert the kwargs dict to named args.
    if 'kwargs' in init_locals:
        kwargs = init_locals.pop('kwargs')
    else:
        kwargs = {}

    args = {
        **init_locals,
        **kwargs
    }
    params = ', '.join([f'{k}={repr(v)}' for k, v in args.items()])
    return f'{class_name}({params})'


def make_params_repr(args):
    '''
    >>> print(make_params_repr({ "uid": 1, "name": "My Dataset"}))
    uid=1, name='My Dataset'
    '''
    params = ', '.join([f'{k}={repr(v)}' for k, v in args.items()])
    return params
