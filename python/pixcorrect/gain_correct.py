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
from pixcorrect import decaminfo

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

        saturate = 0.
        for amp in decaminfo.amps:
            sec = section2slice( image['DATASEC'+amp])
            gain = image['GAIN'+amp]
            image.data[sec]*=gain

            # Adjust the weight or variance image if present:
            if image.weight is not None:
                image.weight[sec] *= 1./(gain*gain)
            if image.variance is not None:
                image.variance[sec] *= (gain*gain)

            # Adjust keywords 
            image['GAIN'+amp] = image['GAIN'+amp] / gain
            image['SATURAT'+amp] = image['SATURAT'+amp] * gain
            saturate = max(saturate, image['SATURAT'+amp])

        # The SATURATE keyword is assigned to maximum of the two amps.
        image['SATURATE'] = saturate

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
