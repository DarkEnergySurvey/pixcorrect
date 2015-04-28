#!/usr/bin/env python
"""Gain Correct image (convert pixel values from ADU to electrons)
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, items_must_match
from despyfits.DESImage import DESImage, DESImageCStruct, weight_dtype, section2slice
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'addweight'

class AddWeight(PixCorrectImStep):
    description = "Add and populate a weight plane to an image"
    step_name = config_section

    @classmethod
    def __call__(cls, image, dome):
        """Add a weight plane 

        :Parameters:
            - `image`: the DESImage for weight plane to be added 

        Applies "in place"
        """
 
        logger.info('Adding Weight Image')

        if image.weight is None:
            image.init_weight()
            # Check that dome and data are from same CCD
            try:
                items_must_match(image, dome, 'CCDNUM')
            except:
                return 1
            # Transform the sky image into a variance image
            data=image.data
            var = np.array(data, dtype = weight_dtype)
            for amp in decaminfo.amps:
                sec = section2slice(image['DATASEC'+amp])
                invgain = (image['FLATMED'+amp]/image['GAIN'+amp]) / dome.data[sec]
                var[sec] += image['RDNOISE'+amp]**2 * invgain
                var[sec] *= invgain
            # Add noise from the dome flat shot noise, if present
            if dome.weight is not None:
                var += data * data / (dome.weight*dome.data * dome.data)
            elif dome.variance is not None:
                var += data * data * dome.variance / (dome.data * dome.data)

            image.weight = 1.0/var
            logger.info('Finished building a weight plane')
        else:
            logger.info('Weight plane already present... skipping.')

        ret_code=0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for addition of a weight plane.

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        logger.info('Weight will be added to %s' % image)
    
        flat_fname = config.get(cls.step_name, 'flat')
        logger.info('Reading flat correction from %s'% flat_fname)
        flat = DESImage.load(flat_fname)
        ret_code = cls.__call__(image,flat)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """
        parser.add_argument('--flat', default=None,
                            help='Dome flat correction image')

add_weight = AddWeight()

# internal functions & classes

if __name__ == '__main__':
    add_weight.main()
