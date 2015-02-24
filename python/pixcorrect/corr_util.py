"""Common for image correction steps"""

# imports
from functools import wraps
from ConfigParser import SafeConfigParser 
from argparse import ArgumentParser
import logging
from contextlib import closing
from os import path
import time
import platform
import ctypes

import pyfits
import numpy as np


from pixcorrect import proddir

# constants

logger = logging.getLogger('pixcorrect')

# exception classes

class ArrayShapeException(Exception):
    pass

# This exception to be thrown when a call to an external library
# returns a non-zero value, indicating an error of some sort
class LibraryException(Exception):
    def __init__(self, value):
        self.value = value

# interface functions

def load_shlib(shlib_name):
    """Load a shared library
    """

    lib_ext = {'Linux': 'so',
               'Darwin': 'dylib'}
    fname = shlib_name + '.' + lib_ext[platform.system()]
    try:
        shlib = ctypes.CDLL(fname)
    except KeyError:
        raise RuntimeError, ("Unknown platform: " + platform.system())
        
    return shlib
    

# A decorator that uses a FITS keyword to determine whether a step
# has already been performed, and skips it if it has.
def do_once(arg_idx, fits_keyword):
    def make_do_once_wrapper(f):
        @wraps(f)
        def f_wrapper(*args, **kwargs):
            hdu = args[arg_idx]
            if fits_keyword in hdu.header.keys():
                done = hdu.header[fits_keyword]
            else:
                done = False
            if not done:
                result = f(*args, **kwargs)
                hdu.header[fits_keyword]=time.asctime(time.localtime())
            else:
                result = None
                logger.warning("Skipping " + f.__name__ +  " (" + fits_keyword + " already set)")
            return result
        return f_wrapper
    return make_do_once_wrapper

# contract conditions

def match_array_shapes(arglist, func_args=None, func_kwargs=None, func_result=None):
    """A dbc contract condition that forces multiple arguments to be the same shape 

    A function suitable for passing as the first argument to the
    precondition, postcondition, or invariantcondition decorators in
    the dbc module.

    :Parameters:
        -`arglist`: a list of argument indexes to check

    Raises an ImageShapeException if all of the arguments in arglist are not
    numpy arrays with the same shape.

    """
    ref_shape = func_args[arglist[0]].shape
    for i in arglist[1:]:
        if func_args[arglist[i]].shape != ref_shape:
            raise ImageShapeException() 

def no_lib_error(func_args=None, func_kwargs=None, func_result=None):
    """A dbc contract condition that checks that the returned value is 0

    A function suitable for passing as the first argument to the
    precondition, postcondition, or invariantcondition decorators in
    the dbc module.

    Raises a LibraryException if the wrapped function returns a 
    non-zero value. This is useful when wrapping C libraries that
    follow the converntion that returning a non-zero value 
    indicates an error.
    """
    if func_result != 0:
        raise LibraryException(func_result) 

# classes
# internal functions & classes
