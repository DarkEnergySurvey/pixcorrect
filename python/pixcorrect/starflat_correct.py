#!/usr/bin/env python
"""Apply a star flat correction to a raw DES image.
Uses all the code from flat_combine but changes only the FITS
keywords used to register actions
"""

from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once
from despyfits.DESImage import DESImage
from despyfits import maskbits
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo
from pixcorrect.flat_correct import FlatCorrect

# Which section of the config file to read for this step
config_section = 'starflat'

class StarFlatCorrect(FlatCorrect):
    description = "Apply a flat field correction to an image"
    step_name = config_section

    @classmethod
    @do_once(1,'DESSTAR')
    def __call__(cls, image, flat_im):
        """Apply a flat field correction to an image

        :Parameters:
            - `image`: the DESImage to flatten
            - `flat_im`:  the star flat correction image to apply

        """
        ret_code = cls._doit(image, flat_im)

        if flat_im.sourcefile is None:
            image.write_key('STARFIL', 'UNKNOWN', comment='Star flat correction file')
        else:
            image.write_key('STARFIL', path.basename(flat_im.sourcefile), comment='Star flat correction file')

        return ret_code

starflat_correct = StarFlatCorrect()

# internal functions & classes

if __name__ == '__main__':
    starflat_correct.main()
