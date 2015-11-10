#!/usr/bin/env python
"""Do image-by-image pixel level corrections
"""

# imports
from functools import partial
import ctypes
import sys

import numpy as np
import pyfits

from despyfits.DESImage import DESImage, DESBPMImage, section2slice 

from pixcorrect import corr_util
from pixcorrect import imtypes
from pixcorrect.dbc import precondition, postcondition
from pixcorrect.corr_util import logger

from pixcorrect.bias_correct import bias_correct
from pixcorrect.make_mask import make_mask
from pixcorrect.fix_columns import fix_columns
from pixcorrect.linearity_correct import linearity_correct
from pixcorrect.gain_correct import gain_correct
from pixcorrect.flat_correct_cp import flat_correct_cp
from pixcorrect.null_weights import null_weights
from pixcorrect.sky_compress import sky_compress
from pixcorrect.sky_subtract import sky_subtract
from pixcorrect.bf_correct import bf_correct
from pixcorrect.add_weight import add_weight
from pixcorrect.starflat_correct import starflat_correct
from pixcorrect import bfinfo
from pixcorrect import skyinfo
from pixcorrect import decaminfo

from pixcorrect.PixCorrectDriver import PixCorrectMultistep

class PixCorrectCP(PixCorrectMultistep):
    config_section = "pixcorrect_cp"
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


    @classmethod
    def _check_return(cls,retval):
        """
        Exit the program if the retval is nonzero.
        """
        if retval!=0:
            sys.exit(retval)
        return
    
    def __call__(self):
        """Do image-by-image pixel level corrections
        """
        # All the code here, asside from one call for each step, should 
        # be assiciated with shoveling data between steps. Everything else should
        # take inside the code for its respective step.

        # Get the science image
        self.sci = DESImage.load(self.config.get('pixcorrect_cp','in'))

        # Bias subtraction
        if self.do_step('bias'):
            self._check_return(bias_correct(self.sci, self.bias))
        self.clean_im('bias')

        # Linearization
        if self.do_step('lincor'):
            lincor_fname=self.config.get('pixcorrect_cp','lincor')
            self._check_return(linearity_correct(self.sci,lincor_fname))

        # Make the mask plane and mark saturated pixels.  Note that flags
        # are set to mark saturated pixels and keep any previously existing mask bits.
        if self.do_step('bpm'):
            self._check_return(make_mask(self.sci,
                                         self.bpm,
                                         saturate=True,
                                         clear=False))

        flat_gaincorrect = self.config.getboolean('pixcorrect_cp','flat_gaincorrect')
        # get gains ahead of time so that jump can be removed from a flat with no gain correction
        gain_preserve={}
        if (flat_gaincorrect):
            tmp_gains={}
            avg_gain=0.0
            for amp in decaminfo.amps:
                tmp_gains[amp] = self.sci['GAIN'+amp]
                avg_gain=avg_gain+tmp_gains[amp]
            for amp in decaminfo.amps:
                gain_preserve[amp]=2.0*tmp_gains[amp]/avg_gain
#            print avg_gain
#            print gain_preserve

        if self.do_step('gain'):
            self._check_return(gain_correct(self.sci))
            
        # B/F correction
        if self.do_step('bf'):
            bf_fname = self.config.get('pixcorrect_cp', 'bf')
            self._check_return(bf_correct(self.sci,
                                          bf_fname,
                                          bfinfo.DEFAULT_BFMASK))
            
        # If done with the BPM; let python reclaim the memory
        if not self.do_step('fixcol'):
            self.clean_im('bpm')


        # Flat field
        if self.do_step('flat'):
#            allow_mismatch = self.config.get('pixcorrect_cp','flat_gaincorrect')
            print "flat_gaincorrect: ",flat_gaincorrect
#            for amp in decaminfo.amps:
#                self.flat[gain_preserve[amp]['sec']]*=gain_preserve[amp]['cor']
            self._check_return(flat_correct_cp(self.sci, self.flat, gain_preserve))
            if not self.do_step('sky'):
                self.clean_im('flat')

        # Fix columns
        if self.do_step('fixcols'):
            self._check_return(fix_columns(self.sci, self.bpm))
            self.clean_im('bpm')

        # Make mini-sky image
        if self.do_step('mini'):
            mini = self.config.get('pixcorrect_cp','mini')
            blocksize = self.config.getint('pixcorrect_cp','blocksize')
            self._check_return(sky_compress(self.sci,
                                            mini,
                                            blocksize,
                                            skyinfo.DEFAULT_SKYMASK))

        # Subtract sky and make weight plane - forcing option to do "sky-only" weight
        if self.do_step('sky'):
            sky_fname = self.config.get('pixcorrect_cp','sky')
            fit_fname = self.config.get('pixcorrect_cp','skyfit')
            self._check_return(sky_subtract(self.sci,
                                            fit_fname,
                                            sky_fname,
                                            'sky',
                                            self.flat))
            if not self.do_step('addweight'):
                self.clean_im('flat')
        
        # Star flatten
        if self.do_step('starflat'):
            self._check_return(starflat_correct(self.sci, self.starflat))
        self.clean_im('starflat')

        ### Do add_weight before null_weight step, else it will overwrite the nulls
        if self.do_step('addweight'):
            self._check_return(add_weight(self.sci, self.flat))
        self.clean_im('flat')

        # This new call should take care of both --resaturate and --null_mask
        if self.do_step('null_mask') or self.do_step('resaturate'):
            # We need to fix the step_name if we want to call 'step_run'
            null_weights.__class__.step_name = self.config_section
            logger.info("Running null_weights")
            self._check_return(null_weights.step_run(self.sci,self.config))

        out_fname = self.config.get('pixcorrect_cp', 'out')
        self.sci.save(out_fname)

        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to pixcorrect_cp driver
        """
        parser.add_argument('--bias', default=None,
                            help='Bias correction image')
        parser.add_argument('--lincor', default=None, 
                            help='linearity correction Table')
        parser.add_argument('--bf', default=None, 
                            help='brighter/fatter correction Table')
        parser.add_argument('--gain', action='store_true', default=False,
                            help='convert ADU to e- using gain values in hdr')
        parser.add_argument('--bpm', default=None, 
                            help='bad pixel mask filename')
        parser.add_argument('--flat', default=None,
                            help='Dome flat correction image')
        parser.add_argument('--flat_gaincorrect', action='store_true', default=False,
                            help='Flag to gain correction of an uncorrected flat')
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
        parser.add_argument('--addweight', action='store_true', default=False,
                            help='Add a weight map to the image if none exists')

        # Adds --resaturate and --null_mask from null_weights class
        null_weights.add_step_args(parser)

        return

if __name__ == '__main__':
    PixCorrectCP.main()
