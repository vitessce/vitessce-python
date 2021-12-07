
def make_repr(init_locals, class_name=None):
    '''
    >>> from .wrappers import MultiImageWrapper
    >>> orig = MultiImageWrapper('IMAGE_WRAPPERS', foo='bar')
    >>> orig_repr = repr(orig)
    >>> print(orig_repr)
    MultiImageWrapper(image_wrappers='IMAGE_WRAPPERS', use_physical_size_scaling=False, foo='bar')
    >>> evalled = eval(orig_repr)
    >>> assert orig_repr == repr(evalled)
    '''
    if 'self' in init_locals:
        del init_locals['self']
    if '__class__' in init_locals:
        clazz = init_locals.pop('__class__')  # Requires superclass to be initialized.
        if class_name is None:
            class_name = clazz.__name__
    
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
    params = ', '.join([f'{k}={repr(v)}' for k, v in args.items()])
    return params