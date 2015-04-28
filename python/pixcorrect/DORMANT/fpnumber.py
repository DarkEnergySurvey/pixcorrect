#!/usr/bin/env python
"""Mask saturated pixels
"""

# imports
import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESFocalPlaneImages import DESFocalPlaneImages
from despyfits.DESFocalPlaneImages import FocalPlaneCStructArray
from pixcorrect.PixCorrectDriver import PixCorrectFPStep

# constants

# Which section of the config file to read for this step
config_section = 'fpnumber'

# exception classes
# interface functions
# classes

# Lowest level access to the C library function
fpnumber_lib = load_shlib('libfpnumber')
fpnumber_c = fpnumber_lib.fpnumber
fpnumber_c.restype = ctypes.c_int
fpnumber_c.argtypes = [ctypes.POINTER(FocalPlaneCStructArray)]

class FPNumber(PixCorrectFPStep):
    description = "Fill in data with HDU numbers"
    step_name = config_section

    @classmethod
    def __call__(cls, images):
        """Replace data values with HDU numbers
        
        :Parameters:
            - `images`: the array of images

        Applies the correction "in place"
        """
        logger.info('Filling in data with HDU numbers')
        c_call_status = fpnumber_c(images.cstruct)
        logger.info('Finished')
        return c_call_status

fpnumber = FPNumber()

# internal functions & classes

if __name__ == '__main__':
    fpnumber.main()
