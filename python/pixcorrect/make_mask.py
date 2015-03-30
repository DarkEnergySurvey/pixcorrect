#!/usr/bin/env python
"""Apply BPM to mask plane and/or flag saturated pixels
"""

from os import path
import numpy as np
import time
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESImage, DESBPMImage, section2slice
from despyfits.maskbits import *
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'mask'

class MakeMask(PixCorrectImStep):
    description = "Build mask plane, setting appropriate bits from BPM and/or saturated pixels"
    step_name = config_section

    DEFAULT_SATURATE = False
    DEFAULT_CLEAR = False
    
    @classmethod
    def __call__(cls, image, bpm_im, saturate, clear):
        """Create or update the mask plane of an image

        :Parameters:
            - `image`: the DESImage to operate upon.  Mask plane is created if absent
            - `bpm_im`: the DESBPMImage with the bad pixel mask. Skips BPM step if None
            - `saturate`: boolean flag indicating whether to set BADPIX_SATURATE flags
            - `clear`: if True, clear pre-existing mask.  If False, or new bits with old.

        """

        if image.mask is None:
            image.init_mask()
        elif clear:
            image.mask.fill(0)

        ret_code = 0
        if bpm_im is None and saturate is False:
            logger.warning('Null operation requested in make_mask')
            return ret_code

        if bpm_im is not None:
            # Check for header keyword of whether it's been done
            kw = 'DESBPM'
            if kw in image.header.keys():
                logger.warning('Skipping BPM application ('+kw+' already set)')
            else:
                logger.info('Applying BPM')
                # Mark the unusable data
                bitmask = BPMDEF_FLAT_MIN | \
                    BPMDEF_FLAT_MAX | \
                    BPMDEF_FLAT_MASK | \
                    BPMDEF_BIAS_HOT | \
                    BPMDEF_BIAS_WARM | \
                    BPMDEF_BIAS_MASK | \
                    BPMDEF_BIAS_COL | \
                    BPMDEF_CORR | \
                    BPMDEF_WACKY_PIX
                # ??? Add FUNKY_COL to the list of unusable pixels?
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_BPM

                # Mark edge pixels and bad amplifier with their own bits
                bitmask = BPMDEF_EDGE
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_EDGE
                bitmask = BPMDEF_BADAMP
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_BADAMP

                # Mark slightly dodgy pixels
                bitmask = BPMDEF_FUNKY_COL | \
                    BPMDEF_TAPE_BUMP
                # ??? Is this what to do with FUNKY_COL ???
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_SUSPECT
              
                image[kw] = time.asctime(time.localtime())
                image.write_key(kw, time.asctime(time.localtime()),
                                comment = 'Construct mask from BPM')
                if bpm_im.sourcefile is None:
                    image.write_key('BPMFIL', 'UNKNOWN', comment='BPM file used to build mask')
                else:
                    image.write_key('BPMFIL', path.basename(bpm_im.sourcefile), comment='BPM file used to build mask')
                        
                logger.debug('Finished applying BPM')

        if saturate:
            # Check for header keyword of whether it's been done
            kw = 'DESSAT'
            if kw in image.header.keys():
                logger.warning('Skipping saturation check ('+kw+' already set)')
            else:
                logger.info('Flagging saturated pixels')
                nsat = 0
                for amp in decaminfo.amps:
                    sec = section2slice(image['DATASEC'+amp])
                    sat = image['SATURAT'+amp]
                    satpix = image.data[sec]>=sat
                    image.mask[sec][satpix] |= BADPIX_SATURATE
                    nsat += np.count_nonzero(satpix)

                image.write_key(kw, time.asctime(time.localtime()),
                                comment = 'Flag saturated pixels')
                image.write_key('NSATPIX',nsat,
                                comment='Number of saturated pixels')
                
                logger.debug('Finished flagging saturated pixels')

        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the BPM

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        if config.has_option(cls.step_name, 'bpm'):
            bpm_fname = config.get(cls.step_name, 'bpm')
            logger.info('reading BPM from %s' % bpm_fname)
            bpm_im = DESBPMImage.load(bpm_fname)
        else:
            bpm_im = None

        if config.has_option(cls.step_name, 'saturate'):
            saturate = config.getboolean(cls.step_name, 'saturate')
        else:
            saturate = DEFAULT_SATURATE

        if config.has_option(cls.step_name, 'clear'):
            clear = config.getboolean(cls.step_name, 'clear')
        else:
            clear = DEFAULT_CLEAR
            
        ret_code = cls.__call__(image, bpm_im, saturate, clear)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('-b', '--bpm', 
                            help='bad pixel mask filename (optional)')
        parser.add_argument('--saturate', action='store_true',
                            help='Flag saturated pixels')
        parser.add_argument('--clear', action='store_true',
                            help='Clear any pre-existing mask bits')

make_mask = MakeMask()

# internal functions & classes

if __name__ == '__main__':
    make_mask.main()
