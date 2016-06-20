#!/usr/bin/env python
"""
"""

from os import path
import sys
import numpy as np
import time
import fitsio

from pixcorrect.corr_util import logger
from despyfits import maskbits
import logging
from argparse import ArgumentParser

class CoaddPrepare:
    description = "Prepare a coadd tile image for SExtractor by interpolating vertically across" \
        "pixels with zero weight, giving a nonzero weight, and marking them in a mask plane."

    DEFAULT_MINCOLS = 1   # Narrowest feature to interpolate
    DEFAULT_MAXCOLS = None  # Widest feature to interpolate.  None means no limit.
    DEFAULT_WEIGHT_THRESHOLD = 0. # Upper limit for weight on invalid pixel
    DEFAULT_WEIGHT_VALUE = -1. # Value placed into masked pixels; negative will copy neighbors' weight
    MASK_VALUE = 1 # Value placed into mask plane for invalid pixel

    @classmethod
    def main(cls):
        description = "Prepare coadd file for SExtractor by interpolating across low-weight pixels," \
          " and adding a mask plane that flags the interpolated pixels"
    
        # Get arguments
        parser = ArgumentParser(description=description)
        parser.add_argument('-l', '--log', 
                        default="", 
                        help="the name of the logfile")
        parser.add_argument('-v', '--verbose', action="count", 
                        help="be verbose")
        parser.add_argument('-i', '--infile', 
                        default=None,
                        help='input coadd image file name')
        parser.add_argument('-o', '--outfile', 
                        default=None,
                        help='output coadd image file name')
        parser.add_argument('--min_cols',
                        default=None, 
                        help='minimum height of region to interpolate')
        parser.add_argument('--max_cols',
                        default=None, 
                        help='maximum height of region to interpolate')
        parser.add_argument('--weight_threshold',
                        default=CoaddPrepare.DEFAULT_WEIGHT_THRESHOLD, 
                        help='Maximum weight value to interpolate over')
        parser.add_argument('--weight_value',
                        default=CoaddPrepare.DEFAULT_WEIGHT_VALUE, 
                        help='Weight value assigned to interpolated pixels, <0 to use neighbors')
    
        args = parser.parse_args()
    
        # Set up logger
        if args.log is not None and len(args.log)>0:
            logging.basicConfig(filename=args.log,
                            format="%(asctime)s %(levelname)s:\t%(message)s",
                            level=logging.WARNING)
            sh = logging.StreamHandler()
            sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s:\t%(message)s"))
            logger.addHandler(sh)
        else:
            logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                            level=logging.WARNING)

        logger = logging.getLogger()
        if args.verbose > 0:
            verbosity = logging.INFO if args.verbose==1 else logging.DEBUG
            logger.setLevel(verbosity)

        # Call routine
        ret_val = coadd_prepare(args.infile, args.outfile,
                            min_cols=args.min_cols,
                            max_cols=args.max_cols,
                            weight_threshold = float(args.weight_threshold),
                            weight_value = float(args.weight_value))
        sys.exit(ret_val)
        
    @classmethod
    def __call__(cls, imageIn, imageOut,
                 min_cols=DEFAULT_MINCOLS,
                 max_cols=DEFAULT_MAXCOLS,
                 weight_threshold = DEFAULT_WEIGHT_THRESHOLD,
                 weight_value = DEFAULT_WEIGHT_VALUE):
        """
        Add a mask plane to imageIn and set bits wherever the weight value is <= given
        threshold.  Then interpolate the data plane along columns to replace masked pixels.
        Set weight to weight_value at the masked pixels too.  Masked pixels along edges are
        not interpolated, and are left with input weight.

        :Parameters:
            - `imageIn`: filename of input image
            - `imageOut`: filename of output image
            - `min_cols`: Minimum width of region to be interpolated.
            - `max_cols`: Maximum width of region to be interpolated.
            - `weight_threshold`: Upper bound for weight values to mark as bad
            - `weight_value`: New weight value for bad pixels, enter <0 to use neighbor weights
        """
 
        logger.info('Preparing coadd {:s} for SExtractor'.format(imageIn))

        # Read weight plane and science plane
        sci,scihdr = fitsio.read(imageIn, ext=0, header=True)
        wt,wthdr = fitsio.read(imageIn, ext=1, header=True)

        # Make mask plane
        mask = wt <= float(weight_threshold)

        # Identify column runs to interpolate, start by marking beginnings of runs
        work = np.array(mask)
        work[1:,:] = np.logical_and(mask[1:,:], ~mask[:-1,:])
        xstart,ystart = np.where(work.T)

        # Now ends of runs
        work = np.array(mask)
        work[:-1,:] = np.logical_and(mask[:-1,:], ~mask[1:,:])
        xend, yend = np.where(work.T)
        yend = yend + 1   # Make the value one-past-end

        # If we've done this correctly, every run has a start and an end, on same col
        if not np.all(xstart==xend):
            logger.error("Logic problem, xstart and xend not equal.")
            print xstart,xend ###
            return 1

        # Narrow our list to runs of the desired length range and
        # not touching the edges
        use = yend-ystart >= min_cols
        if max_cols is not None:
            use = np.logical_and(yend-ystart<=max_cols, use)
        use = np.logical_and(ystart>0, use)
        use = np.logical_and(yend<mask.shape[0], use)
        ystart = ystart[use]
        yend = yend[use]
        xstart = xstart[use]

        # Assign mean of top and bottom to runs, and fill in weight plane
        for run in range(len(xstart)):
            sci[ystart[run]:yend[run],xstart[run]] = \
              0.5*(sci[ystart[run]-1,xstart[run]] +
                   sci[yend[run],xstart[run]])
            if weight_value<0:
                fill_weight = 0.5*(wt[ystart[run]-1,xstart[run]] \
                  + wt[yend[run],xstart[run]])
            else:
                fill_weight = weight_value
            wt[ystart[run]:yend[run],xstart[run]] = fill_weight

        # Add to image history
        scihdr['HISTORY'] =time.asctime(time.localtime()) + \
            ' coadd_prepare with weight threshold {:f}'.format(weight_threshold)
                      
        # Write out all three planes
        mask = np.array(mask, dtype=np.int16)*cls.MASK_VALUE
        logger.debug('Writing output images')
        with fitsio.FITS(imageOut, mode=fitsio.READWRITE, clobber=True) as ff:
            ff.write(sci, extname='SCI', header=scihdr, clobber=True)
            ff.write(mask, extname='MSK')
            ff.write(wt, extname='WGT', header=wthdr)
        
        logger.debug('Finished coadd_prepare')

        ret_code=0
        return ret_code


coadd_prepare = CoaddPrepare()

# Call from command line:

if __name__ == '__main__':
    CoaddPrepare.main()
