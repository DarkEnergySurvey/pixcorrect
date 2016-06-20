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

# Which section of the config file to read for this step
config_section = 'rowinterp'

class RowInterp(PixCorrectImStep):
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

    @classmethod
    def __call__(cls, image,
                 min_cols=DEFAULT_MINCOLS,
                 max_cols=DEFAULT_MAXCOLS,
                 interp_mask=DEFAULT_INTERP_MASK,
                 invalid_mask=DEFAULT_INVALID_MASK):
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
 
        logger.info('Interpolating along rows')

        if image.mask is None:
            logger.error('Input image does not have mask')
            return 1

        interpolate = np.array(image.mask & interp_mask, dtype=bool)

        # Make arrays noting where a run of bad pixels starts or ends
        # Then make arrays has_?? which says whether left side is valid
        # and an array with the value just to the left/right of the run.
        work = np.array(interpolate)
        work[:,1:] = np.logical_and(interpolate[:,1:], ~interpolate[:,:-1])
        ystart,xstart = np.where(work)

        work = np.array(interpolate)
        work[:,:-1] = np.logical_and(interpolate[:,:-1], ~interpolate[:,1:])
        yend, xend = np.where(work)
        xend = xend + 1   # Make the value one-past-end

        # If we've done this correctly, every run has a start and an end.
        if not np.all(ystart==yend):
            logger.error("Logic problem, ystart and yend not equal.")
            return 1

        # Narrow our list to runs of the desired length range
        use = xend-xstart >= min_cols
        if max_cols is not None:
            use = np.logical_and(xend-xstart<=max_cols, use)
        xstart = xstart[use]
        xend = xend[use]
        ystart = ystart[use]

        # Now determine which runs have valid data at left/right
        xleft = np.maximum(0, xstart-1)
        has_left = ~np.array(image.mask[ystart,xleft] & invalid_mask, dtype=bool)
        has_left = np.logical_and(xstart>=1,has_left)
        left_value = image.data[ystart,xleft]

        xright = np.minimum(work.shape[1]-1, xend)
        has_right = ~np.array(image.mask[ystart,xright] & invalid_mask, dtype=bool)
        has_right = np.logical_and(xend<work.shape[1],has_right)
        right_value = image.data[ystart,xright]
        
        # Assign right-side value to runs having just right data
        for run in np.where(np.logical_and(~has_left,has_right))[0]:
            image.data[ystart[run],xstart[run]:xend[run]] = right_value[run]
            image.mask[ystart[run],xstart[run]:xend[run]] |= maskbits.BADPIX_INTERP
        # Assign left-side value to runs having just left data
        for run in np.where(np.logical_and(has_left,~has_right))[0]:
            image.data[ystart[run],xstart[run]:xend[run]] = left_value[run]
            image.mask[ystart[run],xstart[run]:xend[run]] |= maskbits.BADPIX_INTERP

        # Assign mean of left and right to runs having both sides
        for run in np.where(np.logical_and(has_left,has_right))[0]:
            image.data[ystart[run],xstart[run]:xend[run]] = \
              0.5*(left_value[run]+right_value[run])
            image.mask[ystart[run],xstart[run]:xend[run]] |= maskbits.BADPIX_INTERP

        # Add to image history
        image['HISTORY'] =time.asctime(time.localtime()) + \
            ' row_interp over mask 0x{:04X}'.format(interp_mask)
                      
        logger.debug('Finished interpolating')

        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for row interpolation

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        if config.has_option(cls.step_name, 'min_cols'):
            min_cols = config.getint(cls.step_name, 'min_cols')
        else:
            min_cols = cls.DEFAULT_MINCOLS
        if config.has_option(cls.step_name, 'max_cols'):
            max_cols = config.getint(cls.step_name, 'max_cols')
        else:
            max_cols = cls.DEFAULT_MAXCOLS
        if config.has_option(cls.step_name, 'interp_mask'):
            interp_mask = maskbits.parse_badpix_mask(config.get(cls.step_name,
                                                                'interp_mask'))
        else:
            interp_mask = maskbits.parse_badpix_mask(cls.DEFAULT_INTERP_MASK)
        if config.has_option(cls.step_name, 'invalid_mask'):
            invalid_mask = maskbits.parse_badpix_mask(config.get(cls.step_name,
                                                                 'invalid_mask'))
        else:
            invalid_mask = maskbits.parse_badpix_mask(cls.DEFAULT_INVALID_MASK)

        ret_code = cls.__call__(image, min_cols, max_cols, interp_mask, invalid_mask)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--min_cols', nargs=1, 
                            default=None, 
                            help='minimum width of region to interpolate')
        parser.add_argument('--max_cols', nargs=1, 
                            default=None, 
                            help='maximum width of region to interpolate')
        parser.add_argument('--interp_mask', nargs=1, 
                            default=None, 
                            help='bitmask for MSK plane defining pixels to interpolate')
        parser.add_argument('--invalid_mask', nargs=1, 
                            default=None, 
                            help='bitmask for MSK plane defining pixels'
                            ' unusable for interpolation')
        return

row_interp = RowInterp()

# internal functions & classes

if __name__ == '__main__':
    row_interp.main()
