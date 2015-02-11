#!/usr/bin/env python
"""Apply a bad pixel map (pbm) to a DES image 
"""

# imports
import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESBPMImage, DESImageCStruct
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# constants

# Which section of the config file to read for this step
config_section = 'fixcol'

# exception classes
# interface functions
# classes

# Lowest level access to the C library function
fixcol_lib = load_shlib('libfixcol')
fix_cols_c = fixcol_lib.fixCol
fix_cols_c.restype = ctypes.c_int
fix_cols_c.argtypes = [DESImageCStruct, DESImageCStruct]

class FixCols(PixCorrectImStep):
    description = "Fix cols marked in a bad pixel mask"
    step_name = config_section

    @classmethod
    def __call__(cls, image, bpm_im):
        """Apply a bad pixel mask to a DES image

        :Parameters:
            - `image`: the DESImage to which the BPM is to be applied
            - `bpm_im`: the DESImage with the bad pixel mask

        Applies the correction "in place"
        """

        logger.info('Fixing columns marked by BPM')
        ret_code = fix_cols_c(bpm_im.cstruct, image.cstruct)
        logger.debug('Finished fixing columns')
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for masking columns based on the mask

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

fix_cols = FixCols()

# internal functions & classes

if __name__ == '__main__':
    fix_cols.main()
