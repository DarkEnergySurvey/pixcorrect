#!/usr/bin/env python
"""Mask saturated pixels
"""

# imports
import ctypes
from os import path
import numpy as np
from pixcorrect.dbc import postcondition
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from pixcorrect.corr_util import do_once, no_lib_error
from despyfits.DESImage import DESImage, DESImageCStruct
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# constants

# Which section of the config file to read for this step
config_section = 'mask_saturation'

# exception classes
# interface functions
# classes

# Lowest level access to the C library function
masksatr_lib = load_shlib('libmasksatr')
mask_saturation_c = masksatr_lib.mask_saturation
mask_saturation_c.restype = ctypes.c_int
mask_saturation_c.argtypes = [DESImageCStruct, ctypes.POINTER(ctypes.c_int)]

class MaskSaturation(PixCorrectImStep):
    description = "Mark saturated pixels in the mask"
    step_name = config_section

    @classmethod
    @do_once(1,'DESSAT')
    @postcondition(no_lib_error)
    def __call__(cls, image):
        """Mark saturated pixels in the mask of an image
        
        :Parameters:
            - `image`: the DESImage in which to mask saturated pixels

        Applies the correction "in place"
        """
        logger.info('Masking saturated pixels')
        num_saturated = ctypes.c_int()
        c_call_status = mask_saturation_c(image.cstruct, num_saturated)
        logger.info('Masked %d pixels as saturated' % num_saturated.value)
        return c_call_status

mask_saturation = MaskSaturation()

# internal functions & classes

if __name__ == '__main__':
    mask_saturation.main()
