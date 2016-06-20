#!/usr/bin/env python
"""Do image-by-image pixel level corrections
"""

# imports
from functools import partial
import ctypes
import sys

import numpy as np
import pyfits

from despyfits.DESFocalPlaneImages import DESFocalPlaneImages

from pixcorrect import corr_util
from pixcorrect import imtypes
from pixcorrect.dbc import precondition, postcondition
from pixcorrect.corr_util import logger

from pixcorrect.nullop_fp import nullop_fp
from pixcorrect.apply_bpm import apply_bpm
from pixcorrect.override_bpm import override_bpm
from pixcorrect.fix_cols import fix_cols
from pixcorrect.mask_saturation import mask_saturation
from pixcorrect.PixCorrectDriver import PixCorrectMultistep



class PixCorrectFP(PixCorrectMultistep):
    config_section = "pixcorrect_fp"
    step_name = config_section
    description = 'Do image-by-image pixel level corrections'

    def image_data(self, image_name):
        """Return a DESFocalPlaneImages object for the configured images

        :Parameters:
            -`image_name`: the type of image to return

        @returns: the object of class DESFocalPlaneImage
        """
        # If we already have the data, return it
        if image_name in self._image_data:
            im = self._image_data[image_name]
        else:
            # If we don't already have the data, load it
            fname = self.config.get(self.config_section, image_name)
            im = DESFocalPlaneImages.load(fname)
            logger.info('Reading %s image from %s' % (image_name, fname))
            self._image_data[image_name] = im

        return im

    def __call__(self):
        """Do image-by-image pixel level corrections
        """
        # All the code here, asside from one call for each step, should 
        # be assiciated with shoveling data between steps. Everything else should
        if self.do_step('nullop_fp'):
            nullop_fp(self.sci)

        out_fname_template = self.config.get(self.config_section, 'out')
        self.sci.save(out_fname_template)

        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to pixcorrect driver
        """
        parser.add_argument('--nullop_fp', action='store_true',
                                      help='perform a nullop')

if __name__ == '__main__':
    PixCorrectFP.main()
