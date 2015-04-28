#!/usr/bin/env python
"""Override a bad pixel map (pbm) to a DES image 
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct
from pixcorrect.apply_bpm import ApplyBPM

# Which section of the config file to read for this step
config_section = 'obpm'

# Lowest level access to the C library function
bpm_lib = load_shlib('libbpm')
obpm_c = bpm_lib.obpm
obpm_c.restype = ctypes.c_int
obpm_c.argtypes = [DESImageCStruct, DESImageCStruct]

class OverrideBPM(ApplyBPM):

    description = "Override a bad pixel mask"
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
        ret_code = obpm_c(image.cstruct, bpm_im.cstruct)
        logger.debug('Finished applying BPM')
        return ret_code

override_bpm = OverrideBPM()

# internal functions & classes

if __name__ == '__main__':
    override_bpm.main()
