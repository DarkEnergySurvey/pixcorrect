#!/usr/bin/env python
"""Apply BPM to mask plane and/or flag saturated pixels
"""

from os import path
import numpy as np
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESImage
from despyfits.maskbits import *
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'nullweight'

class NullWeightsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class NullWeights(PixCorrectImStep):
    description = "Set weights to zero based on bitmask, and put high values into saturated pixels"
    step_name = config_section

    DEFAULT_RESATURATE = False
    DEFAULT_BITMASK = '0'
    
    @classmethod
    def __call__(cls, image, bitmask, resaturate):
        """Create or update the mask plane of an image

        :Parameters:
            - `image`: the DESImage to operate upon.  Mask plane is created if absent
            - `bitmask`: Integer or list of BADPIX bit names that, when set in mask image,
                         will put weight=0 for that pixel.
            - `resaturate`: if True, set data for every pixel with BADPIX_SATURATE set
                          to a value above the SATURATE keyword

        """

        if image.mask is None:
            raise NullWeightsError('Mask is missing in image')

        if bitmask!=0:
            logger.info('Nulling weight image from mask bits')
            
            if image.weight is None and image.variance is None:
                raise NullWeightsError('Weight is missing in image')
            weight = image.get_weight()
            kill = np.array( image.mask & bitmask, dtype=bool)
            weight[kill] = 0.
            image.header.write_history(time.asctime(time.localtime()) + \
                                       ' Null weights with mask 0x{:04X}'.format(bitmask))
            logger.debug('Finished nulling weight image')
            
        if resaturate:
            logger.info('Re-saturating pixels from mask bits')
            sat = np.array( image.mask & BADPIX_SATURATE, dtype=bool)
            image.data[sat] = 1.01 * image['SATURATE']
            image.header.write_history(time.asctime(time.localtime()) + \
                                       ' Set saturated pixels to {:.0f}')
            logger.debug('Finished nulling weight image')


        ret_code = 0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the BPM

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        if config.has_option(cls.step_name, 'bitmask'):
            bitmask = parse_badpix_mask(config.get(cls.step_name, 'bitmask'))
        else:
            bitmask = parse_badpix_mask(cls.DEFAULT_BITMASK)

        if config.has_option(cls.step_name, 'resaturate'):
            resaturate = config.getboolean(cls.step_name, 'resaturate')
        else:
            resaturate = cls.DEFAULT_RESATURATE

        ret_code = cls.__call__(image, bitmask, resaturate)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('--resaturate', action='store_true',
                            help='Put saturated value in BADPIX_SATURATE pixels')
        parser.add_argument('--bitmask', default=cls.DEFAULT_BITMASK,
                            help='Clear any pre-existing mask bits')

null_weights = NullWeights()

# internal functions & classes

if __name__ == '__main__':
    null_weights.main()
