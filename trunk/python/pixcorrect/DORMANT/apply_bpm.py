#!/usr/bin/env python
"""Apply a bad pixel map (pbm) to a DES image 
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESBPMImage, DESImageCStruct
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'bpm'

# Lowest level access to the C library function
bpm_lib = load_shlib('libbpm')
bpm_c = bpm_lib.bpm
bpm_c.restype = ctypes.c_int
bpm_c.argtypes = [DESImageCStruct, DESImageCStruct]

class ApplyBPM(PixCorrectImStep):
    description = "Apply a bad pixel mask"
    step_name = config_section

    @classmethod
    def __call__(cls, image, bpm_im):
        """Apply a bad pixel mask to a DES image

        :Parameters:
            - `image`: the DESImage to which the BPM is to be applied
            - `bpm_im`: the DESImage with the bad pixel mask

        Applies the correction "in place"
        """

        logger.info('Applying BPM')
        if image.mask is None:
            image.init_mask()
        ret_code = bpm_c(image.cstruct, bpm_im.cstruct)
        logger.debug('Finished applying BPM')
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the BPM

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        bpm_fname = config.get(cls.step_name, 'bpm')
        logger.info('reading BPM from %s' % bpm_fname)
        bpm_im = DESBPMImage.load(bpm_fname)
    
        ret_code = cls.__call__(image, bpm_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('-b', '--bpm', nargs=1, 
                            default=None, 
                            help='bad pixel mask filename')

apply_bpm = ApplyBPM()

# internal functions & classes

if __name__ == '__main__':
    apply_bpm.main()
