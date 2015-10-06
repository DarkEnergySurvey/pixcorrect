#!/usr/bin/env python
"""Apply BPM to mask plane and/or flag saturated pixels
"""

from os import path
import numpy as np
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESImage
from despyfits import maskbits
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo
import time

# Which section of the config file to read for this step
config_section = 'nullweight'

class NullWeightsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class NullWeights(PixCorrectImStep):
    description = "Set weights to zero based on null_mask, and put high values into saturated pixels"
    step_name = config_section

    DEFAULT_RESATURATE = False
    DEFAULT_NULL_MASK = '0'
    
    @classmethod
    def __call__(cls, image, null_mask, resaturate):
        """Create or update the mask plane of an image

        :Parameters:
            - `image`: the DESImage to operate upon.  Mask plane is created if absent
            - `null_mask`: Integer or list of BADPIX bit names that, when set in mask image,
                         will put weight=0 for that pixel.
            - `resaturate`: if True, set data for every pixel with BADPIX_SATURATE set
                          to a value above the SATURATE keyword

        """

        if image.mask is None:
            raise NullWeightsError('Mask is missing in image')

        if null_mask!=0:
            logger.info('Nulling weight image from mask bits')
            
            if image.weight is None and image.variance is None:
                raise NullWeightsError('Weight is missing in image')
            weight = image.get_weight()
            kill = np.array( image.mask & null_mask, dtype=bool)
            weight[kill] = 0.
            image['HISTORY'] =time.asctime(time.localtime()) + \
                              ' Null weights with mask 0x{:04X}'.format(null_mask)
            logger.debug('Finished nulling weight image')
            
        if resaturate:
            logger.info('Re-saturating pixels from mask bits')
            sat = np.array( image.mask & maskbits.BADPIX_SATURATE, dtype=bool)
            try:
                saturation_level = image['SATURATE']
            except (ValueError,KeyError):
                # If there is no SATURATE, try taking max of amps
                maxsat = 0.
                try:
                    for amp in decaminfo.amps:
                        maxsat = max(maxsat, image['SATURAT'+amp])
                except:
                    logger.error('SATURATx header keywords not found')
                    raise NullWeightsError('SATURATx header keywords not found')
                saturation_level = maxsat
                logger.warning('Taking SATURATE as max of single-amp SATURATx values')
                
            image.data[sat] = 1.01 * saturation_level
            image['HISTORY'] = time.asctime(time.localtime()) + \
                              ' Set saturated pixels to {:.0f}'.format(saturation_level)
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
        # Debug config
        #print "###### CONFIG"
        #for section in config.sections():
        #    print section
        #    for option in config.options(section):
        #        print " ", option, "=", config.get(section, option)
        #print "########"
        
        if config.has_option(cls.step_name, 'null_mask'):
            null_mask = maskbits.parse_badpix_mask(config.get(cls.step_name, 'null_mask'))
        else:
            null_mask = maskbits.parse_badpix_mask(cls.DEFAULT_NULL_MASK)

        if config.has_option(cls.step_name, 'resaturate'):
            resaturate = config.getboolean(cls.step_name, 'resaturate')
        else:
            resaturate = cls.DEFAULT_RESATURATE

        ret_code = cls.__call__(image, null_mask, resaturate)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('--resaturate', action='store_true',
                            help='Put saturated value in BADPIX_SATURATE pixels')
        parser.add_argument('--null_mask', default=cls.DEFAULT_NULL_MASK,
                            help='Names of mask bits to null (or an integer mask)')

null_weights = NullWeights()

# internal functions & classes

if __name__ == '__main__':
    null_weights.main()
