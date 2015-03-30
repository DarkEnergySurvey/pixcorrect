#!/usr/bin/env python
"""Apply a linearity correction to a DES image 
"""

import ctypes
from os import path
import numpy as np
import fitsio
from scipy import interpolate
from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once
from despyfits.DESImage import DESImage, section2slice
from despyfits.DESFITSInventory import DESFITSInventory
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'lincor'

class LinearityCorrect(PixCorrectImStep):
    description = "Apply a linearity correction"
    step_name = config_section

    @classmethod
    @do_once(1,'DESLINC')
    def __call__(cls, image, fname_lincor):
        """Apply a linearity correction

        :Parameters:
            - `image`: the DESImage to determine and apply an ovescan correction
            - `fname_lincor`: the linearity correction FITS table (contains look-up tables)

        Applies the correction "in place"
        """

#
#       Discover the HDU in the linearity correction FITS table that contains data for a specific CCD
#
        fits_inventory = DESFITSInventory(fname_lincor)
        lincor_hdu=fits_inventory.ccd_hdus(image['CCDNUM'])
        if (len(lincor_hdu) != 1):
            if (len(lincor_hdu) == 0):
                logger.error('Unable to locate HDU in %s containing linearity correction for CCDNUM %d. Aborting!'.format(fname_lincor,image['CCDNUM']))
            else:
                logger.error('Found multiple HDUs in %s containing linearity correction for CCDNUM %d. Aborting!'.format(fname_lincor,image['CCDNUM']))
            raise

        logger.info('Reading Linearity Correction from %s' % (fname_lincor))
        cat_fits=fitsio.FITS(fname_lincor,'r')
        cat_hdu=lincor_hdu[0]
        cols_retrieve=["ADU","ADU_LINEAR_A","ADU_LINEAR_B"]
        CAT=cat_fits[cat_hdu].read(columns=cols_retrieve)
#
#        If columns do not get put into CAT in a predefined order then these utilities
#        may be needed.  RAG has them and can implement... left this way for now since it 
#        currently duplicates imcorrect exactly
#
#        CATcol=cat_fits[cat_hdu].get_colnames()
#        cdict=MkCatDict(CATcol,cols_retrieve)

#
#       Define the correction being made.
#
        nonlinear=[]
        linearA=[]
        linearB=[]
        for row in CAT:
            nonlinear.append(row[0])
            linearA.append(row[1])
            linearB.append(row[2])
        nonlinear=np.array(nonlinear)
        linearA=np.array(linearA)
        linearB=np.array(linearB)
        interpA = interpolate.interp1d(nonlinear, linearA, kind='linear', copy=True)
        interpB = interpolate.interp1d(nonlinear, linearB, kind='linear', copy=True)
        logger.info('Applying Linearity Correction')

#
#       Slice over the datasecs for each amplifier.
#       Apply the correction
#
        seca = section2slice( image['DATASECA'])
        secb = section2slice( image['DATASECB'])

        # Only fix pixels that are in the range of the nonlinearity table
        in_range=np.logical_and(image.data[seca]>=np.min(nonlinear),
                                image.data[seca]<=np.max(nonlinear))
        image.data[seca][in_range]=interpA(image.data[seca][in_range])

        in_range=np.logical_and(image.data[secb]>=np.min(nonlinear),
                                image.data[secb]<=np.max(nonlinear))
        image.data[secb][in_range]=interpB(image.data[secb][in_range])

        image.write_key('LINCFIL', path.basename(fname_lincor), comment='Nonlinearity correction file')
        
        ret_code=0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the BPM

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        fname_lincor = config.get(cls.step_name, 'lincor')
        logger.info('Linearity correction will be applied to %s' % image)
    
        ret_code = cls.__call__(image, fname_lincor)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('--lincor', nargs=1, type=str, default=None, 
                            help='Linearity Correction Table')


linearity_correct = LinearityCorrect()

# internal functions & classes

if __name__ == '__main__':
    linearity_correct.main()
