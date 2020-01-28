#!/usr/bin/env python3

# $Id: cti.py 47952 2019-01-03 21:04:53Z rgruendl $
# $Rev:: 47952                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2019-01-03 15:04:53 #$:  # Date of last commit.

"""Perform CTI Search and Mask on image
"""

from pixcorrect import cti_utils as cti
from pixcorrect import lightbulb_utils as lb
#from pixcorrect.corr_util import logger, do_once
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectImStep
#from despyfits.DESImage import DESImage, DESImageCStruct, section2slice, data_dtype

# Which section of the config file to read for this step
config_section = 'cticheck'

class CTIcheck(PixCorrectImStep):
    description = "Search for indication that a known CTI (Charge Transfer Inefficiency) is active and mask"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """
        This is currently written with instance of CTI (Charge Transfer Inefficiency) in mind (that occuring
        for CCD=41 during DES Y6).  It may be generalizable if further cases occur (but not all the parts have
        yet been written with a generalized case in mind).  When CTI is detected the entire amplifier in
        question will be masked with BADPIX_BADAMP
        """
#       A simple dictionary with parameters for the only known case of CTI
#           Currently explist is set to encompass DES Y6 (20180815 and beyond (expnum>765533)
#           This could be tightened to a range as no CTI has been detected after November 2018 but
#           it has not yet been systematicall watched for.
        CTI = {41: {'amp': 'B', 'explist': '765533-'}}

        if image['CCDNUM'] in CTI:

#
#           Currently can re-use the function developed for lightbulb checking
#
            check_for_light = lb.check_lightbulb_explist(image['EXPNUM'], CTI[image['CCDNUM']]['explist'])
            if check_for_light:
                logger.info(' CTI: Expnum={:d}, CCDNUM={:d}, in proscribed range checking for CTI'.format(
                    image['EXPNUM'], image['CCDNUM']))
                ctiDict = cti.check_cti(image, CTI[image['CCDNUM']], verbose=1)
#
#               Current criterion:
#                   Looks for horizontal striping in image (with large deficits in counts that are not
#                       associated with an edge-bleed.
#                   Examines auto-correlation for lags in the x-direction at 5, 7, and 15 pixel offsets
#                       and compares to lags obtained from measurments in the diaganol direction.
#                   Looks for evidence of excessive power in the ratio between x-direction and diagnol sets
#                       that indicative that charge is bleeding in the x-direction.
#
                if ctiDict['isCTI']:
                    image = cti.mask_cti(image, CTI[image['CCDNUM']], ctiDict, verbose=1)
                    logger.info(' CTI: Detected CTI for Exp={:d}, CCD={:d}, Amp={:s}'.format(image['EXPNUM'], image['CCDNUM'], CTI[image['CCDNUM']]['amp']))
                    image.write_key('DES_CTI', 'Masked DATASEC{:s}'.format(CTI[image['CCDNUM']]['amp']))

        logger.debug('Finished checking and applying mask CTI')
        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for check and masking of CTI

        :Parameters:
            - `image`: the DESImage on which to operate
#            - `config`: the configuration from which to get other parameters

        """
        logger.info('CTI check %s' % image)

        ret_code = cls.__call__(image)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """

cticheck = CTIcheck()

# internal functions & classes

if __name__ == '__main__':
    cticheck.main()
