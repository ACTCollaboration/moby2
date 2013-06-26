from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import warnings
import functools
import inspect

class MobyDeprecationWarning(Warning):
    pass

def deprecated(func):
    """
    This is a decorator which can be used to mark functions as
    deprecated. It will result in a warning being emitted when the
    function is used.  Adapted from the PythonDecoratorLibrary.
    """
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
            "call to deprecated function {}.".format(func.__name__),
            category=MobyDeprecationWarning,
            filename=func.__code__.co_filename,
            lineno=func.__code__.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func

def deprecatedmethod(func):
    """
    This is a decorator which can be used to mark class member
    functions as deprecated.
    """
    @functools.wraps(func)
    def new_func(self, *args, **kwargs):
        call_loc = inspect.stack()[-1]
        warnings.warn_explicit(
            "call to deprecated class member function {}.{}".format(
                str(self.__class__.__name__),func.__name__),
            category=MobyDeprecationWarning,
            filename=call_loc[1],
            lineno=call_loc[2]
        )
        return func(self, *args, **kwargs)
    return new_func
