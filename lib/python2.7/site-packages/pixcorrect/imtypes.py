"""Tools to verify image types"""

# imports
from pixcorrect.dbc import make_args_check

# constants

sci_shape = (2048, 4096)

# exception classes

class ImageTypeException(Exception):
    pass

class ImageWrongHeader(ImageTypeException):
    def __init__(self, kw, read_value, expected_value):
        self.kw = kw
        self.read_value = read_value
        self.expected_value = expected_value

class ImageWrongShape(ImageTypeException):
    def __init__(self, read_shape, expected_shape):
        self.read_shape = read_shape
        self.expected_shape = expected_shape

# interface functions

# I could have just made one class, and made the checkers for the different
# image types instances of this class. However, using different classes and
# staticmetheds makes it easier to allow arbitrary checking code to be added
# for specific image types, while this would be awkward if they were all
# instances of the same class.

class ImageTypeChecker(object):
    shape = None
    kwdict = {}
    keywords = []

    @classmethod
    def check(cls, hdu):
        for kw in cls.kwdict:
            if not kw in hdu.header:
                raise ImageWrongHeader(kw, 'missing', None)
            if hdu.header[kw] != cls.kwdict[kw]:
                raise ImageWrongHeader(kw, cls.kwdict[kw], hdu.header[kw])

        for kw in cls.keywords:
            if not kw in hdu.header:
                raise ImageWrongHeader(kw, 'missing', None)

        if cls.shape is not None:
            if hdu.data.shape != cls.shape:
                raise ImageWrongShape(hdu.data.shape, cls.shape)

class Type0ImageChecker(ImageTypeChecker):
    shape = sci_shape

class Type1ImageChecker(Type0ImageChecker):
    kwdict = {'FORALL': True,
              'TREE': 'oak'}

class Type3ImageChecker(Type0ImageChecker):
    keywords = ['BARCOEFF']    

#
# Make design-by-contract (dbc) decorators from our checking classes
#

args_type0 = make_args_check(Type0ImageChecker.check)
args_type1 = make_args_check(Type1ImageChecker.check)
args_type3 = make_args_check(Type3ImageChecker.check)

# internal functions & classes
