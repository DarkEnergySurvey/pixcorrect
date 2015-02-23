#!/usr/bin/env python
"""Apply a bias correction to a raw DES image 
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'bias'

# Lowest level access to the C library function
biascorrect_lib = load_shlib('libbiascorrect')
bias_c = biascorrect_lib.bias_c
bias_c.restype = ctypes.c_int
bias_c.argtypes = [DESImageCStruct, DESImageCStruct]

class BiasCorrect(PixCorrectImStep):
    description = "Apply a bias correction to an image"
    step_name = config_section

    @classmethod
    def __call__(cls, image, bias_im):
        """Apply a bias correction to an image

        :Parameters:
            - `image`: the DESImage to apply a bias correction
            - `bias_im`:  the bias correction image to apply

        Applies the correction "in place"
        """
 
        logger.info('Applying Bias')
        ret_code = bias_c(image.cstruct, bias_im.cstruct)
        logger.debug('Finished applying Bias')
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `bias`: the bias image to apply

        """

        bias_fname = config.get(cls.step_name, 'bias')
        logger.info('reading Bias from %s'% bias_fname)
        bias_im = DESImage.load(bias_fname)
    
        ret_code = cls.__call__(image, bias_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the bias correction
        """
        parser.add_argument('--bias', nargs=1, default=None,
                            help='Bias correction image')

bias_correct = BiasCorrect()

# internal functions & classes

if __name__ == '__main__':
    bias_correct.main()
