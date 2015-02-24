#!/usr/bin/env python
"""Apply a flat correction to a raw DES image 
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'flat'

# Lowest level access to the C library function
flatcorrect_lib = load_shlib('libflatcorrect')
flat_c = flatcorrect_lib.flat_c
flat_c.restype = ctypes.c_int
flat_c.argtypes = [DESImageCStruct, DESImageCStruct]

class FlatCorrect(PixCorrectImStep):
    description = "Apply a flat field correction to an image"
    step_name = config_section

    @classmethod
    def __call__(cls, image, flat_im):
        """Apply a flat field correction to an image

        :Parameters:
            - `image`: the DESImage to apply a bias correction
            - `flat_im`:  the flat correction image to apply

        Applies the correction "in place"
        """
 
        logger.info('Applying Flat')
        ret_code = flat_c(image.cstruct, flat_im.cstruct)
        logger.debug('Finished applying Flat')
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `flat`: the bias image to apply

        """

        flat_fname = config.get(cls.step_name, 'flat')
        logger.info('Reading flat correction from %s'% flat_fname)
        flat_im = DESImage.load(flat_fname)
    
        ret_code = cls.__call__(image, flat_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--flat', nargs=1, default=None,
                            help='Flat field correction image')

flat_correct = FlatCorrect()

# internal functions & classes

if __name__ == '__main__':
    flat_correct.main()
