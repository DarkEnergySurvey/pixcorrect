#!/usr/bin/env python
"""
Implement John Marriner's correctable-column fixer.
I have made a few changes to the C fixCols() function:
* Using clipped mean instead of median to estimate sky levels.
* Slightly different way to choose comparison columns.  Shouldn't matter
* Set maximum distance that a comparison column can be.
* Instead of refusing to fix any column that has any saturated pixel in it, I exclude
  from the correction any pixels that are saturated.
* When calculating comparison columns, restrict to same rows that have target pixels in them.
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once, items_must_match
from despyfits.DESImage import DESDataImage, DESImage, DESBPMImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect.clippedMean import clippedMean
from despyfits import maskbits

# Which section of the config file to read for this step
config_section = 'fixcolumns'

class FixColumnsError(Exception):
    """
    Error class for problems in fixing columns
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class FixColumns(PixCorrectImStep):
    description = "Fix the correctable columns"
    step_name = config_section
    
    CORR = maskbits.BPMDEF_CORR  # BPM flag for correctable pixels
    REJECT = maskbits.BADPIX_BPM + \
      maskbits.BADPIX_SATURATE +\
      maskbits.BADPIX_BADAMP + \
      maskbits.BADPIX_SUSPECT  # bitmask for comparison pixels
    MINIMUM_PIXELS = 100  # Smallest number of pixels in column to correct
    CLIP_SIGMA = 4    # Rejection threshold for mean statistics

    @classmethod
    def _valid_pix(cls, image, bpm, icol):
        """
        Return boolean array saying which pixels in column icol are useful for sky stats
        """
        use = (bpm.mask[:,icol] & cls.REJECT)==0
        use &= ~np.isinf(image.data[:,icol])
        use &= ~np.isnan(image.data[:,icol])
        return use
    
    @classmethod
    @do_once(1,'DESFIXC')
    def __call__(cls, image, bpm):
        """
        Find and fix correctable columns in the image as indicated by the BPMDEF_CORR
        bit being set in the bpm image.  The algorithm is taken from John Marriner's
        fixCol() in old mask_utils.c.  The affected pixels in the column have a constant
        added to them that makes their median value equal that in neighboring columns.

        :Parameters:
            - `image`: DESImage to fix.
            - `bpm`: DESBPMImage for this CCD
        """
 
        logger.info('Fixing columns')

        NEIGHBORS = 10  # Number of comparison columns to seek
        RANGE = 20  # Farthest away to look for comparison columns
        # Largest allowable fractional difference in variance between the fixable column
        # and its neighbors:
        VAR_TOLERANCE = 0.5
        # John Marriner would not apply correction unless it was this much larger
        # than the statistical noise in the correction:
        MINIMUM_SIGNIFICANCE = 5 
        
        
        if image.mask is None:
            raise FixColumnsError('Input image does not have mask')
        # Check that dome and data are from same CCD
        try:
            items_must_match(image, bpm, 'CCDNUM')
        except:
            return 1

        # A "fixable" column will have CORR flag set at either start or end of column
        fixable = np.where(np.logical_or(bpm.mask[0,:] & cls.CORR,
                                         bpm.mask[-1,:] & cls.CORR))[0]

        for icol in fixable:
            # Which pixels in the column are fixable?
            # They need to have only the CORR flag set, and be finite, and not saturated.
            coldata = image.data[:,icol]
            colbpm = bpm.mask[:,icol]
            ignore = np.logical_or( colbpm & ~cls.CORR, np.isinf(coldata))
            ignore |= np.isnan(coldata)
            ignore |= image.mask[:,icol] & maskbits.BADPIX_SATURATE
            use_rows = np.logical_and(colbpm & cls.CORR, ~ignore)
            if np.count_nonzero(use_rows) < cls.MINIMUM_PIXELS:
                logger.info("Not enough pixels to fix column {:d}".format(icol))
                continue

            # Get a robust estimate of mean level in target column
            col_mean, col_var, col_n = clippedMean(coldata[use_rows], cls.CLIP_SIGMA)

            # Now want to collect stats on up to NEIGHBORS nearby columns
            norm_stats = []
            ilow = icol
            ihigh = icol
            low_limit = max(icol - RANGE,0)
            high_limit = min(icol + RANGE, image.data.shape[1]-1)
            while len(norm_stats) < NEIGHBORS and (ilow>low_limit or ihigh<high_limit):
                while ilow>low_limit:
                    # get stats from next useful column to left:
                    ilow-=1
                    if ilow in fixable:
                        continue
                    use = cls._valid_pix(image, bpm, ilow)
                    use &= use_rows
                    if np.count_nonzero(use) < cls.MINIMUM_PIXELS:
                        continue
                    norm_stats.append(clippedMean(image.data[:,ilow][use],cls.CLIP_SIGMA))
                    break
                while ihigh<high_limit:
                    # get stats from next useful column to right:
                    ihigh+=1
                    if ihigh in fixable:
                        continue
                    use = cls._valid_pix(image, bpm, ihigh)
                    use &= use_rows
                    if np.count_nonzero(use) < cls.MINIMUM_PIXELS:
                        continue
                    norm_stats.append(clippedMean(image.data[:,ihigh][use],cls.CLIP_SIGMA))
                    break
            if len(norm_stats) < NEIGHBORS:
                # Don't fix the column if we did not get comparison columns
                logger.info('Not enough comparison columns to fix col {:d}'.format(icol))
                continue

            # Calculate the weighted mean estimate of mean neighbor cols
            mean = np.array([i[0] for i in norm_stats])
            var = np.array([i[1] for i in norm_stats])
            wt = np.array([i[2] for i in norm_stats]) / var

            # Do not apply correction if the target column's variance is much
            # different from the comparison columns
            norm_var = np.sum(var*wt)/np.sum(wt)
            if np.abs(col_var - norm_var) > VAR_TOLERANCE * norm_var:
                logger.info('Too much variance to fix col {:d}'.format(icol))
                continue

            correction = np.sum(mean*wt)/np.sum(wt) - col_mean
            correction_var = 1./np.sum(wt) + col_var/col_n
            # Marriner does not apply correction if it's insignificant:
            if correction*correction < correction_var * MINIMUM_SIGNIFICANCE**2:
                print 'correction:::',correction ##
                logger.info('Insignificant correction for ' \
                             'column {:d} by {:f}'.format(icol,float(correction)))
                continue

            # Apply correction:
            image.data[:,icol][use_rows] += correction
            # Promote the corrected pixels from useless to just imperfect:
            image.mask[:,icol][use_rows] &= ~maskbits.BADPIX_BPM
            image.mask[:,icol][use_rows] |= maskbits.BADPIX_SUSPECT
            print 'correction:::',correction ##
            logger.info('Corrected column {:d} by {:f}'.format(icol,float(correction)))

        if bpm.sourcefile is None:
            image.write_key('FIXCFIL', 'UNKNOWN', comment='BPM file for fixing columns')
        else:
            image.write_key('FIXCFIL', path.basename(bpm.sourcefile), comment='BPM file for fixing columns')

        logger.debug('Finished fixing columns')

        

        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for sky subtraction

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        bpm_fname = config.get(cls.step_name, 'bpm')
        logger.info('reading BPM from %s' % bpm_fname)
        bpm_im = DESBPMImage.load(bpm_fname)
    
        ret_code = cls.__call__(image, bpm_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('-b', '--bpm', nargs=1, 
                            default=None, 
                            help='bad pixel mask filename')
        return

fix_columns = FixColumns()

# internal functions & classes

if __name__ == '__main__':
    fix_columns.main()
