#!/usr/bin/env python
"""Gain Correct image (convert pixel values from ADU to electrons)
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, section2slice, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'gain'

class GainCorrect(PixCorrectImStep):
    description = "Convert image units from ADU to Electrons"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """Convert pixel values from ADU to electrons, including weight or variance
        image and critical keywords.

        :Parameters:
            - `image`: the DESImage for pixel values to be converted

        Applies the correction "in place"
        """
 
        logger.info('Gain Correcting Image')

        seca = section2slice( image['DATASECA'])
        secb = section2slice( image['DATASECB'])

        gaina = image['GAINA']
        gainb = image['GAINB']
        
        image.data[seca]*=gaina
        image.data[secb]*=gainb

        if image.weight is not None:
            image.weight[seca] *= 1./(gaina*gaina)
            image.weight[secb] *= 1./(gainb*gainb)
        if hasattr(image,'variance') and image.variance is not None:
            # The hasattr() check is because DESImage does not yet have variance
            image.variance[seca] *= (gaina*gaina)
            image.variance[secb] *= (gainb*gainb)

        # Adjust keywords 
        image['GAINA'] = image['GAINA'] / gaina
        image['RDNOISEA'] = image['RDNOISEA'] * (gaina*gaina)
        image['SATURATA'] = image['SATURATA'] * gaina
        image['GAINB'] = image['GAINB'] / gainb
        image['RDNOISEB'] = image['RDNOISEB'] * (gainb*gainb)
        image['SATURATB'] = image['SATURATB'] * gainb
        # The SATURATE keyword is assigned to maximum of the two amps.
        # Is this really what we want??? I would think minimum.
        image['SATURATE'] = max(image['SATURATA'],image['SATURATB'])
        
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
