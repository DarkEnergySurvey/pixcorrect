#!/usr/bin/env python
"""Apply a flat correction to a raw DES image 
"""

import ctypes
import sys
import os
#from os import path
import fitsio
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'scaleflat'

# Lowest level access to the C library function
#flatcorrect_lib = load_shlib('libflatcorrect')
#flat_c = flatcorrect_lib.flat_c
#flat_c.restype = ctypes.c_int
#flat_c.argtypes = [DESImageCStruct, DESImageCStruct]

class ScaleFlat(PixCorrectImStep):
    description = "Normalize a single flat field images by a predetermined value"
    step_name = config_section

    @classmethod
    def __call__(cls, image, normfactor, ampborder ):
        """Apply a flat field correction to an image

        :Parameters:
            - `image`: input image (image to be normalized)
            - `normval`: normalization value to use (with respect to SCALMEAN keyword)
            - `ampborder`: area around periphery of AmpRegion to use when calculating statistics

        Applies the correction to each input and writes a separate output file.
        """
 
        logger.info('Normalizing Flat Image')
        scalmean=image['SCALMEAN']
        nfactor=scalmean/normfactor
        nfactor2=nfactor*nfactor
        logger.info('SCALMEAN=%.2f NORMFACTOR=%.2f NORMALIZATION=%.5f' % (scalmean,normfactor,nfactor) )
        image['NORMFACT']=normfactor
#
        image.data*=nfactor            
        if image.weight is not None:
            image.weight*=nfactor2
        elif image.variance is not None:
            image.variance /= nfactor2
#
#       Create keywords that reflect the median value of the flat on each amp.
#
        for amp in decaminfo.amps:
            datasecn=scan_fits_section(image,'DATASEC'+amp)
            datasecn[0]=datasecn[0]+ampborder
            datasecn[1]=datasecn[1]-ampborder
            datasecn[2]=datasecn[2]+ampborder
            datasecn[3]=datasecn[3]-ampborder
            image['FLATMED'+amp]=np.median(image.data[datasecn[2]:datasecn[3]+1,datasecn[0]:datasecn[1]+1])
            
        logger.debug('Finished applying normalization to Flat')
        ret_code=0

        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `normfactorfile`: file containing normalization factor (relative to SCALMEAN) to be used
            - `ampborder`:  the length (in pixels) to exclude around the periphery of each AMP to ignore when calculating statistics

        """

        normfactorfile = config.get(cls.step_name, 'normfactorfile')
        try:
            f_normfactorfile=open(normfactorfile,'r')
        except:
            raise IOError("File not found.  Missing normfactorfile %s " % normfactorfile)
        for line in f_normfactorfile:
            line=line.strip()
            columns=line.split()
            normfactor=float(columns[0])
        f_normfactorfile.close()
        ampborder = config.getint(cls.step_name, 'ampborder')

        ret_code = cls.__call__(image, normfactor, ampborder)

        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--normfactorfile', nargs=1, default=None,
                            help='File containing the normalization factor (relative to SCALMEAN) to be used.')
        parser.add_argument('--ampborder', type=int, default=50,
                            help='Length in pixels around periphery of each amp to ignore when calculating statistics (default=50)')

scale_flat = ScaleFlat()

# internal functions & classes

if __name__ == '__main__':
    scale_flat.main()
