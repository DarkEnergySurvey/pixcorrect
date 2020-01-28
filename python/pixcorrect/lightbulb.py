#!/usr/bin/env python3

# $Id: lightbulb.py 46993 2018-05-10 20:03:29Z rgruendl $
# $Rev:: 46993                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2018-05-10 15:03:29 #$:  # Date of last commit.

"""Lightbulb Search and Mask on image
"""

import numpy as np
from pixcorrect import lightbulb_utils as lb
#from pixcorrect.corr_util import logger, do_once
from pixcorrect.corr_util import logger
#from despyfits.DESImage import DESImage, DESImageCStruct, section2slice, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'lightbulb'

class LightBulb(PixCorrectImStep):
    description = "Search for indication that a known lightbulb defect is active and mask"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """
        This is currently written with a single lightbulb in mind.  It may be generalizable if further cases
        occur (but not all the parts have yet been written with a generalized case in mind).
        Three components are generated depending  on the strength of the lightbulb:
            1) circular mask centered on the bulb with a radius sufficient to reach 1-sigma of the noise level
            2) columnar mask that extends toward the read registers (triggered when lightbulb has central brightness > XX
            3) columnar mask that extends away from the read registers  (triggered when light bulb is saturated)
        """
#       A simple dictionary with parameters for the only known lightbulb
#           Currently explist is set to encompass 20170901 and beyond (expnum>674105)
#           This could be tightened to 693244 (20171101) or even to 694699 (the earliest weak lighbulb found so far)
        LBD = {46: {'explist': '674105-',
                    'xc': 795,
                    'yc': 2620,
                    'rad':500}}

        if image['CCDNUM'] in LBD:
            check_for_light = lb.check_lightbulb_explist(image['EXPNUM'], LBD[image['CCDNUM']]['explist'])
            if check_for_light:
                logger.info('Image CCDNUM=46, in proscribed range checking for lightbulb')
                bulbDict = lb.check_lightbulb(image, LBD[image['CCDNUM']], verbose=1)
#
#               Current criterion:
#               1) saturate pixels detected, median brightness >100,000 and width above the lower end of range
#               2) non-saturated:
#                   requires successful fit (failed fits set g_wid, g_widerr, g_amp, g_amperr to -1)
#                   width (FWHM) in range 70-140 and uncertainty less than 35
#                   amplitude > 0 and uncertainty positive and less than sqrt(amp)
#
                isBulb = False
                bulbDict['bulb'] = 'F'
                if bulbDict['num_sat'] > 0:
                    if bulbDict['g_wid'] > 70. and bulbDict['g_amp'] > 100000:
                        bulbDict['bulb'] = 'S'
                        isBulb = True
                if not isBulb:
                    if bulbDict['g_wid'] >= 70. and bulbDict['g_wid'] <= 140. and \
                       bulbDict['g_widerr'] > 0. and bulbDict['g_widerr'] < 35.:
                        if bulbDict['g_amp'] >= 0.:
                            if bulbDict['g_amperr'] > 0. and bulbDict['g_amperr'] < np.sqrt(bulbDict['g_amp']):
                                bulbDict['bulb'] = 'T'
                                isBulb = True
#
#               If found use masking utility
#
                if isBulb:
                    logger.info(' LIGHTBULB: detected with central brightness: {:.1f}'.format(bulbDict['bulb_sb']))
                    image = lb.mask_lightbulb(image, LBD[image['CCDNUM']], bulbDict, verbose=1)
                    image.write_key('DESBULB', 'Peak SB {:.1f}, Radius {:.1f} '.format(bulbDict['bulb_sb'],
                                                                                       bulbDict['g_wid']),
                                    comment='')

        logger.debug('Finished checking and applying mask for light bulb')
        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for check and masking of light bulb

        :Parameters:
            - `image`: the DESImage on which to operate
#            - `config`: the configuration from which to get other parameters

        """
        logger.info('Light bulb check %s', image)

        ret_code = cls.__call__(image)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """

lightbulb = LightBulb()

# internal functions & classes

if __name__ == '__main__':
    lightbulb.main()
