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
        clazz = init_locals["self"].__class__
    elif class_def is not None:
        clazz = class_def
    else:
        raise ValueError("make_repr could not locate the class definition")

    # Remove self from locals.
    if 'self' in init_locals:
        del init_locals['self']

    # Get the class name.
    class_name = clazz.__name__

    # Remove redundant constructor parameters (when the value equals the default value).
    for k, v in inspect.signature(clazz).parameters.items():
        try:
            if k in init_locals and init_locals[k] == v.default:
                del init_locals[k]
        except:
            # Equality comparison may not be implemented for the value object.
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
