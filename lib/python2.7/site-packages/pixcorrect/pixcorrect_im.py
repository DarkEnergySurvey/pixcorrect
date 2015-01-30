#!/usr/bin/env python
"""Do image-by-image pixel level corrections
"""

# imports
from functools import partial
import ctypes
import sys

import numpy as np
import pyfits

from despyfits.DESImage import DESImage

from pixcorrect import corr_util
from pixcorrect import imtypes
from pixcorrect.dbc import precondition, postcondition
from pixcorrect.corr_util import logger

from pixcorrect.nullop import nullop
from pixcorrect.apply_bpm import apply_bpm
from pixcorrect.override_bpm import override_bpm
from pixcorrect.fix_cols import fix_cols
from pixcorrect.mask_saturation import mask_saturation
from pixcorrect.PixCorrectImDriver import PixCorrectImDriver

config_section = "pixcorrect_im"

class PixCorrectIm(PixCorrectImDriver):
    step_name = config_section
    description = 'Do image-by-image pixel level corrections'
    _image_data = {}
    
    def __init__(self, config):
        self.config = config

    @classmethod
    def run(cls, config):
        config.set(config_section, 'sci', 
                   config.get(config_section, 'in'))
        pix_corrector = cls(config)
        ret_value = pix_corrector()
        return ret_value

    def image_data(self, image_name):
        """Return a DESImage object for a configured image

        :Parameters:
            -`image_name`: the type of image to return

        @returns: the object of class DESImage
        """
        # If we already have the data, return it
        if image_name in self._image_data:
            im = self._image_data[image_name]
        else:
            # If we don't already have the data, load it
            fname = self.config.get('pixcorrect_im', image_name)
            im = DESImage.load(fname)
            logger.info('Reading %s image from %s' % (image_name, fname))
            self._image_data[image_name] = im

        return im

    def __getattr__(self, image_name):
        """Create a shortcut to images using object attributes
        """
        return self.image_data(image_name)

    def clean_im(self, image_name):
        """Let python garbage collect the memory used for an image

        :Parameters:
            -`image_name`: the type of image to clean
        """
        if image_name in self._image_data:
            del self._image_data[image_name]
        

    def __call__(self):
        """Do image-by-image pixel level corrections
        """
        # All the code here, asside from one call for each step, should 
        # be assiciated with shoveling data between steps. Everything else should
        # take inside the code for its respective step.

        # Are we configured to do the named step?
        def do_step(step_name):
            if not self.config.has_option('pixcorrect_im', step_name):
                return False

            try:
                # If the parameter is a boolean, interpret is
                # as an on/off switch
                doit = self.config.getboolean('pixcorrect_im', step_name)
                return doit
            except:
                # Otherwise, interpret it as a value associated with
                # the step, and assume we want to perform the step
                return True

        if do_step('nullop'):
            nullop(self.sci)

        if do_step('bpm'):
            if do_step('override_bpm'):
                override_bpm(self.sci, self.bpm)
            else:
                apply_bpm(self.sci, self.bpm)

        if do_step('fixcol'):
            fix_cols(self.sci, self.bpm)

        # We should be done with the BPM; let python reclaim the memory
        self.clean_im('bpm')

        if do_step('mask_saturation'):
            mask_saturation(self.sci)

        out_fname = self.config.get('pixcorrect_im', 'out')
        self.sci.save(out_fname)

        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to pixcorrect driver
        """
        parser.add_argument('--bpm', nargs=1, 
                                      default=None, 
                                      help='bad pixel mask filename')
        parser.add_argument('--fix_cols', action='store_true',
                                      help='fix bad columns')
        parser.add_argument('--mask_saturation', action='store_true',
                                      help='add saturated pixels to the mask')


if __name__ == '__main__':
    PixCorrectIm.main()
