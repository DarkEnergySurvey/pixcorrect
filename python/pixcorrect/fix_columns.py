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
    # BPM flag for correctable pixels
    CORR = maskbits.BPMDEF_CORR
    # These BPM flags define pixels that are not correctable
    BPMBAD = maskbits.BPMDEF_FLAT_MIN | \
    maskbits.BPMDEF_FLAT_MAX | \
    maskbits.BPMDEF_FLAT_MASK | \
    maskbits.BPMDEF_BIAS_HOT | \
    maskbits.BPMDEF_BIAS_WARM | \
    maskbits.BPMDEF_BIAS_MASK | \
    maskbits.BPMDEF_WACKY_PIX
 
    MINIMUM_PIXELS = 100  # Smallest number of pixels in column to correct
    CLIP_SIGMA = 4    # Rejection threshold for mean statistics

    @classmethod
    def _clippedLine(cls,y,z,doslope,nSigma):
        """
        Perform a straight line fit to data, iteratively clipping outliers > nSigma
        It is implicitly assumed that the slope is small.  The variance is computed
        for slope=0 and outliers are initially clipped with that slope.  This avoids
        the problem of outliers biasing the slope excessively.
        If doslope is true, clipping outliers includes the fitted slope, otherwise it
        does not include the slope.  The slope is always calculated and the mean value
        does not depend on the result for the slope (except in so far as it affects
        clipping of outliers).
        """

        iqdSigma = 1.349
        p25 = np.percentile(z, 25.)
        p50 = np.percentile(z, 50.)
        sigma = (p50-p25)/iqdSigma
        err = nSigma*sigma
        lower = p50 - err
        upper = p50 + err
        mask = np.bitwise_or(z<lower,z>upper) 
        #nrej is number of points rejected in the current pass.
        #Set to arbitrary number>0 for first pass
        nrej = 100
        while nrej:
            yp = y[:][~mask]
            zp = z[:][~mask]
            n = np.size(yp)
            if n<cls.MINIMUM_PIXELS:
                return 0.0, 0.0, -1.0, 0
            avey = np.sum(yp)/n
            yp -= avey
            mean = np.sum(zp)/n
            slope = np.sum(yp*zp)/np.sum(yp*yp)
            if doslope:
                res = z - slope*(y-avey) - mean
            else:
                res = z - mean
            rej = np.absolute(res) > err
            nrej = np.sum(rej & ~mask)
            mask |= rej
            
        gres = res[:][~mask]
        var = np.sqrt(np.sum(gres*gres)/n)
        return mean, slope, var, n

    @classmethod
    def _valid_pix(cls, image, bpm, icol):
        """
        Return boolean array saying which pixels in column icol are useful for sky stats
        """
        #Allow NEAREDGE columns for reference
        use = (image.mask[:,icol] & ~maskbits.BADPIX_NEAREDGE)==0
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
        #Modified 7/7/2016
        #Use clipLine to fit column slope
        #Change NEIGHBORS from 10 to 6
        #Change VAR_TOLERANCE from 0.5 to 0.25
        #Remove lower limit on correction
        #Change mask bit usage
        #Correct all correctable pixels, but use only "good" pixels to compute correction

        logger.info('Fixing columns')
        
        NEIGHBORS = 6  # Number of comparison columns to seek
        RANGE = 12  # Farthest away to look for comparison columns
        # Largest allowable fractional difference in variance between the fixable column
        # and its neighbors:
        VAR_TOLERANCE = 0.25
        #We require the number of sky pixels in the target column to be not less than this
        #fraction of the average number of sky pixels in the reference columns
        #If the pixel values of a column vary in a bi-stable way, the high pixels may be
        #interpreted as "objects" and the high variance may not be noticed.
        COUNT_TOL = 0.85
        
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
        #Just an array that gives the ordinal number of each row in the column (for fitting the slope)
        colord = np.arange(4096)
        #Don't use slope in clipLine
        doslope = False

        for icol in fixable:
            # The fixable column is a hot bias pixel type if COL_BIAS is set
            #hotbias = np.logical_or(bpm.mask[0,icol] & maskbits.BPMDEF_BIAS_COL,
            #                             bpm.mask[-1,icol] & maskbits.BPMDEF_BIAS_COL)
            # Which pixels in the column are fixable?
            # They need to have the CORR flag set (other BPM bits specified by BPMOK are allowed)
            # Checking for valid NAN's or INF's should no longer be necessary, but is harmless
            # Also, we do not use any bad pixels.  
            coldata = image.data[:,icol]
            colbpm = bpm.mask[:,icol]
            ignore = np.logical_or(colbpm & cls.BPMBAD, np.isinf(coldata))
            ignore |= np.isnan(coldata)
            corr_rows = np.logical_and(colbpm & cls.CORR, ~ignore)
            ignore |= image.mask[:,icol] & ~maskbits.BADPIX_BPM
            use_rows = np.logical_and(colbpm & cls.CORR, ~ignore)

            if np.count_nonzero(use_rows) < cls.MINIMUM_PIXELS:
                logger.info("Not enough pixels to fix column {:d}".format(icol))
                continue

            # Get a robust estimate of mean level and slope in target column
            y = colord[use_rows]
            z = coldata[use_rows]
            col_mean, col_slope, col_var, col_n = cls._clippedLine(y,z,doslope,cls.CLIP_SIGMA)
            if col_var <= 0.0:
                logger.info("Error in clipped line fit for column {:d}".format(icol))
                continue
                
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
                    y = colord[use]
                    z = image.data[:,ilow][use]
                    ref_col,ref_slope,ref_var, ref_n = cls._clippedLine(y,z,doslope,cls.CLIP_SIGMA)
                    if ref_var<=0.0: continue
                    norm_stats.append([ref_col,ref_slope,ref_var,ref_n])
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
                    y = colord[use]
                    z = image.data[:,ihigh][use]                
                    ref_col,ref_slope,ref_var,ref_n = cls._clippedLine(y,z,doslope,cls.CLIP_SIGMA)
                    if ref_var<=0.0: continue
                    norm_stats.append([ref_col,ref_slope,ref_var,ref_n])
                    break
                
            if len(norm_stats) < NEIGHBORS:
                # Don't fix the column if we did not get comparison columns
                logger.info('Not enough comparison columns to fix col {:d}'.format(icol))
                continue

            # Calculate the weighted mean estimate of mean neighbor cols
            mean = np.array([i[0] for i in norm_stats])
            var = np.array([i[1] for i in norm_stats])
            wt = np.array([i[2] for i in norm_stats]) / var
            slope = np.array([i[1] for i in norm_stats])
            var = np.array([i[2] for i in norm_stats])
            wt = np.array([i[3] for i in norm_stats]) / var
            nc = np.array([i[3] for i in norm_stats])
            # Do not apply correction if the target column's variance is much
            # different from the comparison columns
            norm_var = np.sum(var*wt)/np.sum(wt)
            if np.abs(col_var - norm_var) > VAR_TOLERANCE * norm_var:
                logger.info('Too much variance to fix column {:d}'.format(icol))
                continue
            #Check that number of target column sky pixels is not much less than
            #the average of the reference columns
            norm_n = np.sum(nc)/np.size(nc)
            if col_n < COUNT_TOL*norm_n:
                logger.info('Too few sky pixels to fix column {:d}'.format(icol))
                continue
 
            #Valid correction.  Calculate correction & error estimate
            norm_mean = np.sum(mean*wt)/np.sum(wt)
            correction = norm_mean - col_mean
            correction_var = 1./np.sum(wt) + col_var/col_n
            # Apply correction:
            image.data[:,icol][corr_rows] += correction
            # Promote the corrected pixels from useless to just imperfect:
            image.mask[:,icol][corr_rows] &= ~maskbits.BADPIX_BPM
            image.mask[:,icol][corr_rows] |= maskbits.BADPIX_FIXED
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
