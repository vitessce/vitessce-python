
def make_repr(init_locals, class_name=None, params_only=False):
    '''
    >>> from .wrappers import MultiImageWrapper
    >>> orig = MultiImageWrapper('IMAGE_WRAPPERS', foo='bar')
    >>> orig_repr = repr(orig)
    >>> print(orig_repr)
    MultiImageWrapper(image_wrappers='IMAGE_WRAPPERS', use_physical_size_scaling=False, foo='bar')
    >>> evalled = eval(orig_repr)
    >>> assert orig_repr == repr(evalled)
    '''
    try:
        del init_locals['self']
    except KeyError:
        pass
    try:
        clazz = init_locals.pop('__class__')  # Requires superclass to be initialized.
        class_name = clazz.__name__
    except KeyError:
        pass
    try:
        kwargs = init_locals.pop('kwargs')
    except KeyError:
        kwargs = {}
    args = {
        **init_locals,
        **kwargs
    }
    params = ', '.join([f'{k}={repr(v)}' for k, v in args.items()])
    if params_only:
        return params
    else:
        return f'{class_name}({params})'