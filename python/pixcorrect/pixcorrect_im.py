#!/usr/bin/env python
"""Do image-by-image pixel level corrections
"""

# imports
from functools import partial
import ctypes
import sys

import numpy as np
import pyfits

from despyfits.DESImage import DESImage, DESBPMImage

from pixcorrect import corr_util
from pixcorrect import imtypes
from pixcorrect.dbc import precondition, postcondition
from pixcorrect.corr_util import logger

from pixcorrect.nullop import nullop
from pixcorrect.nullop_im import nullop_im
from pixcorrect.bias_correct import bias_correct
from pixcorrect.apply_bpm import apply_bpm
from pixcorrect.override_bpm import override_bpm
from pixcorrect.fix_cols import fix_cols
from pixcorrect.mask_saturation import mask_saturation
from pixcorrect.linearity_correct import linearity_correct
from pixcorrect.gain_correct import gain_correct
from pixcorrect.flat_correct import flat_correct
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

class PixCorrectIm(PixCorrectMultistep):
    config_section = "pixcorrect_im"
    step_name = config_section
    description = 'Do image-by-image pixel level corrections'
    _image_types = {'bpm': DESBPMImage}
    
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

            # Get the class of the image we are loading
            try:
                image_class = self._image_types[image_name]
            except KeyError:
                image_class = DESImage

            fname = self.config.get(self.config_section, image_name)
            im = image_class.load(fname)
            logger.info('Reading %s image from %s' % (image_name, fname))
            self._image_data[image_name] = im

        return im


    def __call__(self):
        """Do image-by-image pixel level corrections
        """
        # All the code here, asside from one call for each step, should 
        # be assiciated with shoveling data between steps. Everything else should
        # take inside the code for its respective step.

        if self.do_step('nullop'):
            nullop()

        if self.do_step('nullop_im'):
            nullop_im(self.sci)

        if self.do_step('bias'):
            bias_correct(self.sci, self.bias)
        self.clean_im('bias')

        if self.do_step('bpm'):
            if self.do_step('override_bpm'):
                override_bpm(self.sci, self.bpm)
            else:
                apply_bpm(self.sci, self.bpm)

        if self.do_step('fixcol'):
            fix_cols(self.sci, self.bpm)

        # We should be done with the BPM; let python reclaim the memory
        self.clean_im('bpm')

        if self.do_step('mask_saturation'):
            mask_saturation(self.sci)

	if self.do_step('lincor'):
	    lincor_fname=self.config.get('pixcorrect_im','lincor')
	    linearity_correct(self.sci,lincor_fname)

        if self.do_step('gain'):
            gain_correct(self.sci)

        if self.do_step('flat'):
            flat_correct(self.sci, self.flat)
        self.clean_im('flat')

        out_fname = self.config.get('pixcorrect_im', 'out')
        self.sci.save(out_fname)

        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to pixcorrect driver
        """
        parser.add_argument('--bias', nargs=1, default=None,
                            help='Bias correction image')
        parser.add_argument('--bpm', nargs=1, default=None, 
                            help='bad pixel mask filename')
        parser.add_argument('--fix_cols', action='store_true',
                            help='fix bad columns')
 	parser.add_argument('--lincor', nargs=1, default=None, 
                            help='Linearity Correction Table')
        parser.add_argument('--gain', action='store_true',
                            help='convert ADU to e- using gain values in hdr')
        parser.add_argument('--mask_saturation', action='store_true',
                            help='add saturated pixels to the mask')
        parser.add_argument('--flat', nargs=1, default=None,
                            help='Flat field correction image')

if __name__ == '__main__':
    PixCorrectIm.main()
