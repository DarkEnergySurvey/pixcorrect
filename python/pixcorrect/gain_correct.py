#!/usr/bin/env python
"""Gain Correct image (convert pixel values from ADU to electrons)
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'gain'

class GainCorrect(PixCorrectImStep):
    description = "Convert image units from ADU to Electrons"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """Convert pixel values from ADU to electrons

        :Parameters:
            - `image`: the DESImage for pixel values to be converted

        Applies the correction "in place"
        """
 
        logger.info('Gain Correcting Image')

        print image['GAINA']
        print 5.0*image['GAINA']

        ampsecan = [x-1 for x in scan_fits_section(image, 'AMPSECA')]
        ampsecbn = [x-1 for x in scan_fits_section(image, 'AMPSECB')]
        if ((ampsecan[2]>ampsecan[3])or(ampsecbn[2]>ampsecbn[3])):
            print("AMPSEC keyword detected with reversed order")
        print "AMPSECA: ",ampsecan
        print "AMPSECB: ",ampsecbn
#       order x-ranges
        if (ampsecan[0]>ampsecan[1]):
            xtemp=ampsecan[0]
            ampsecan[0]=ampsecan[1]
            ampsecan[1]=xtemp
        if (ampsecbn[0]>ampsecbn[1]):
            xtemp=ampsecbn[0]
            ampsecbn[0]=ampsecbn[1]
            ampsecbn[1]=xtemp
#       order y-ranges
        if (ampsecan[2]>ampsecan[3]):
            ytemp=ampsecan[2]
            ampsecan[2]=ampsecan[3]
            ampsecan[3]=ytemp
        if (ampsecbn[2]>ampsecbn[3]):
            ytemp=ampsecbn[2]
            ampsecbn[2]=ampsecbn[3]
            ampsecbn[3]=ytemp
        print "AMPSECA(ordered): ",ampsecan
        print "AMPSECB(ordered): ",ampsecbn
#
        gaina=np.empty([image['NAXIS2'],image['NAXIS1']/2],dtype=np.float32)
        gainb=np.empty([image['NAXIS2'],image['NAXIS1']/2],dtype=np.float32)
        gaina[:]=image['GAINA']
        gainb[:]=image['GAINB']
        image.data[ampsecan[2]:ampsecan[3]+1,ampsecan[0]:ampsecan[1]+1]*=gaina
        image.data[ampsecbn[2]:ampsecbn[3]+1,ampsecbn[0]:ampsecbn[1]+1]*=gainb
        gaina=[]
        gainb=[]
#
        logger.debug('Finished applying Gain Correction')
        ret_code=0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the gain correction

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        logger.info('Gain correction will be applied to %s' % image)
    
        ret_code = cls.__call__(image)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """

gain_correct = GainCorrect()

# internal functions & classes

if __name__ == '__main__':
    gain_correct.main()
