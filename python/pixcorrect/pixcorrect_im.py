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

from pixcorrect.bias_correct import bias_correct
from pixcorrect.make_mask import make_mask
from pixcorrect.fix_columns import fix_columns
from pixcorrect.linearity_correct import linearity_correct
from pixcorrect.gain_correct import gain_correct
from pixcorrect.flat_correct import flat_correct
from pixcorrect.null_weights import null_weights
from pixcorrect.sky_compress import sky_compress
from pixcorrect.sky_subtract import sky_subtract
from pixcorrect.bf_correct import bf_correct
from pixcorrect import bfinfo
from pixcorrect import skyinfo

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

        # Get the science image
        self.sci = DESImage.load(self.config.get('pixcorrect_im','in'))
        
        # Bias subtraction
        if self.do_step('bias'):
            bias_correct(self.sci, self.bias)
        self.clean_im('bias')

        # Linearization
        if self.do_step('lincor'):
            lincor_fname=self.config.get('pixcorrect_im','lincor')
            linearity_correct(self.sci,lincor_fname)

        # B/F correction
        if self.do_step('bf'):
            bf_fname = self.config.get('pixcorrect_im', 'bf')
            bf_correct(self.sci, bf_fname, bfinfo.DEFAULT_BFMASK)
            
        # Gain correct
        if self.do_step('gain'):
            gain_correct(self.sci)

        # Make the mask plane and mark saturated pixels.  Note that flags
        # are set to mark saturated pixels and clear any previously existing mask.
        if self.do_step('bpm'):
            make_mask(self.sci, self.bpm, saturate=True, clear=True)

        # We should be done with the BPM; let python reclaim the memory
        if not self.do_step('fixcol'):
            self.clean_im('bpm')

        # Flat field
        if self.do_step('flat') and not self.do_step('noflat'):
            flat_correct(self.sci, self.flat)
            if not self.do_step('sky'):
                self.clean_im('flat')

        # Fix columns
        if self.do_step('fixcols'):
            fix_columns(self.sci, self.bpm)
            self.clean_im('bpm')

        # Make mini-sky image
        if self.do_step('mini'):
            blocksize = self.config.getint('pixcorrect_im','blocksize')
            sky_compress(self.sci, mini, blocksize, skyinfo.DEFAULT_SKYMASK)

        # Subtract sky and make weight plane - forcing option to do "sky-only" weight
        if self.do_step('sky'):
            fit_fname = self.config.get('pixcorrect_im','skyfit')
            sky_subtract(self.sci, fit_fname, sky, weight='sky', sci.flat)
        self.clean('flat')
        
        # Star flatten
        if self.do_step('starflat'):
            flat_correct(self.sci, self.starflat)
        self.clean_im('starflat')

        # Re-saturate pixels (??? Can also null weights at selected mask bits here)
        if self.do_step('resaturate'):
            null_weights(self.sci, bitmask=0, resaturate=True)
        

        out_fname = self.config.get('pixcorrect_im', 'out')
        self.sci.save(out_fname)

        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to pixcorrect driver
        """
        parser.add_argument('--bias', default=None,
                            help='Bias correction image')
        parser.add_argument('--lincor', default=None, 
                            help='linearity correction Table')
        parser.add_argument('--bf', default=None, 
                            help='brighter/fatter correction Table')
        parser.add_argument('--gain', action='store_true',
                            help='convert ADU to e- using gain values in hdr')
        parser.add_argument('--bpm', default=None, 
                            help='bad pixel mask filename')
        parser.add_argument('--flat', default=None,
                            help='Dome flat correction image')
        parser.add_argument('--fixcols', action='store_true',
                            help='fix bad columns')
        parser.add_argument('--mini', default=None,
                            help='compressed sky image filename')
        parser.add_argument('--blocksize', default=skyinfo.DEFAULT_BLOCKSIZE,
                            help='blocksize for compressed sky image')
        parser.add_argument('--sky', default=None,
                            help='Template file for sky subtraction and weight creation.' \
                            ' Requires flat and skyfit files be given.')
        parser.add_argument('--skyfit', default=None,
                            help='MiniDECam file holding sky fit coefficients')
        parser.add_argument('--starflat', default=None,
                            help='Star flat correction image')
        parser.add_argument('--resaturate', action='store_true',
                            help='Raise all saturated pixels above SATURATE value')
        return

if __name__ == '__main__':
    PixCorrectIm.main()
