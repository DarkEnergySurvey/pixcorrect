#!/usr/bin/env python
"""
Fill bad pixels with values to the left and/or right on the same row.
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError
import time

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from despyfits import maskbits
from despyastro import zipper_interp as zipp


# Which section of the config file to read for this step
config_section = 'rowzipper'

class ZipperInterp(PixCorrectImStep):
    description = "Interpolate along rows using mean of pixel values to left " \
      "and/or to right of regions of pixels targeted for interpolation."
    step_name = config_section
    
    DEFAULT_MINCOLS = 1   # Narrowest feature to interpolate
    DEFAULT_MAXCOLS = None  # Widest feature to interpolate.  None means no limit.
    DEFAULT_INTERP_MASK = maskbits.BADPIX_BPM + \
      maskbits.BADPIX_SATURATE +\
      maskbits.BADPIX_CRAY +\
      maskbits.BADPIX_STAR +\
      maskbits.BADPIX_TRAIL +\
      maskbits.BADPIX_EDGEBLEED +\
      maskbits.BADPIX_STREAK   # Mask bits that trigger interpolation
    DEFAULT_INVALID_MASK = maskbits.BADPIX_BPM + \
      maskbits.BADPIX_SATURATE +\
      maskbits.BADPIX_BADAMP +\
      maskbits.BADPIX_CRAY +\
      maskbits.BADPIX_STAR +\
      maskbits.BADPIX_TRAIL +\
      maskbits.BADPIX_EDGEBLEED +\
      maskbits.BADPIX_EDGE +\
      maskbits.BADPIX_STREAK   # Mask bits that invalidate a pixel as a source
      # of data to use in interpolation.
    DEFAULT_BLOCK_SIZE = 1
    DEFAULT_ADD_NOISE = False
    DEFAULT_CLOBBER  = False

    @classmethod
    def __call__(cls, image, mask,
                 interp_mask=DEFAULT_INTERP_MASK,
                 BADPIX_INTERP= maskbits.BADPIX_INTERP,
                 min_cols=DEFAULT_MINCOLS,
                 max_cols=DEFAULT_MAXCOLS,
                 invalid_mask=DEFAULT_INVALID_MASK,
                 add_noise=DEFAULT_ADD_NOISE,
                 clobber = DEFAULT_CLOBBER,
                 block_size = DEFAULT_BLOCK_SIZE,
                 logger=logger):
        """
        Interpolate over selected pixels by inserting average of pixels to left and right
        of any bunch of adjacent selected pixels.  If the interpolation region touches an
        edge, or the adjacent pixel has flags marking it as invalid, than the value at
        other border is used for interpolation.  No interpolation is done if both
        boundary pixels are invalid.

        :Parameters:
            - `image`: DESImage to fix.
            - `min_cols`: Minimum width of region to be interpolated.
            - `max_cols`: Maximum width of region to be interpolated.
            - `interp_mask`: Mask bits that will trigger interpolation
            - `invalid_mask`: Mask bits invalidating a pixel as interpolation source.
        """

        # Pass the locals as kwargs
        kwargs = locals()
        image, mask = zipp.zipper_interp_rows(**kwargs)
        return image,mask

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for row interpolation

        :Parameters:
            - `config`: the configuration from which to get other parameters
        """

        min_cols = config.getint(cls.step_name, 'min_cols')
        max_cols = config.getint(cls.step_name, 'max_cols')
        interp_mask = maskbits.parse_badpix_mask(config.get(cls.step_name,'interp_mask'))
        invalid_mask = maskbits.parse_badpix_mask(config.get(cls.step_name,'invalid_mask'))
        add_noise  = config.getboolean(cls.step_name, 'add_noise')
        clobber  = config.getboolean(cls.step_name, 'clobber')
        block_size  = config.getint(cls.step_name, 'block_size')

        kwargs = locals()

        logger.info("Will run row_zipper function with:")
        for key in kwargs.keys():
            logger.info("--%s %s" % (key,kwargs[key]))

        # Now we call the function
        image.data, image.mask = cls.__call__(image.data, image.mask,
                                              interp_mask=interp_mask,
                                              min_cols=min_cols,
                                              max_cols=max_cols,
                                              invalid_mask=invalid_mask,
                                              add_noise=add_noise,
                                              block_size=block_size,
                                              clobber=clobber)
        return 

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--min_cols', nargs=1, default=cls.DEFAULT_MINCOLS, 
                            help='minimum width of region to interpolate')
        parser.add_argument('--max_cols', nargs=1, default=cls.DEFAULT_MAXCOLS, 
                            help='maximum width of region to interpolate')
        parser.add_argument('--interp_mask', nargs=1, default=cls.DEFAULT_INTERP_MASK, 
                            help='bitmask for MSK plane defining pixels to interpolate')
        parser.add_argument('--invalid_mask', nargs=1, default=cls.DEFAULT_INVALID_MASK,
                            help='bitmask for MSK plane defining pixels unusable for interpolation')
        parser.add_argument("--clobber", action='store_true', default=cls.DEFAULT_CLOBBER,
                            help="Clobber output fits file")
        parser.add_argument("--add_noise", action='store_true', default=cls.DEFAULT_ADD_NOISE,
                            help="Add Poisson Noise to the zipper")
        parser.add_argument("--block_size", type=int,default=cls.DEFAULT_BLOCK_SIZE,
                            help="Block size of zipper in x-direction (row)")
        return

row_zipper = ZipperInterp()

# internal functions & classes

if __name__ == '__main__':
    row_zipper.main()
