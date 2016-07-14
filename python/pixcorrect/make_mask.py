#!/usr/bin/env python
"""Apply BPM to mask plane and/or flag saturated pixels
"""

from os import path
import numpy as np
import time
from pixcorrect.corr_util import logger, items_must_match
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

        #Flag saturated pixels first, if requested
        if saturate:
            # Check for header keyword of whether it's been done
            kw = 'DESSAT'
            if kw in image.header.keys() and not clear:
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

        #Now fill in BPM
        if bpm_im is not None:
            # Check for header keyword of whether it's been done
            kw = 'DESBPM'
            if kw in image.header.keys() and not clear:
                logger.warning('Skipping BPM application ('+kw+' already set)')
            else:
                logger.info('Applying BPM')
                try:
                    items_must_match(image, bpm_im, 'CCDNUM')
                except:
                    return 1

                #====Temporary kluge until we get the new BPMS
                #Replace CORR with BIAS_COL
                #bitmask = BPMDEF_CORR
                #mark = (bpm_im.mask & bitmask) != 0
                #bpm_im.mask[mark] |= BPMDEF_BIAS_COL
                # Clear correctable bits from BPM if any are already set
                #bpm_im.mask -= (bpm_im.mask & BPMDEF_CORR)
                #====End kluge

                # Map the following BPM bits to BADPIX_BPM in the image mask
                bitmask = BPMDEF_FLAT_MIN | \
                    BPMDEF_FLAT_MAX | \
                    BPMDEF_FLAT_MASK | \
                    BPMDEF_BIAS_HOT | \
                    BPMDEF_BIAS_WARM | \
                    BPMDEF_BIAS_MASK | \
                    BPMDEF_BIAS_COL | \
                    BPMDEF_FUNKY_COL | \
                    BPMDEF_WACKY_PIX
                # ERICM Removed BPMDEF_CORR and added FUNKY_COL to the above list 
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_BPM

                # Copy BPM edge pixels to image mask
                bitmask = BPMDEF_EDGE
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_EDGE

                # Copy bad amplifier bits to image mask
                bitmask = BPMDEF_BADAMP
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_BADAMP

                # Copy SUSPECT BPM bits to image mask
                bitmask = BPMDEF_SUSPECT
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_SUSPECT

                # Copy NEAREDGE BPM bits to image mask
                bitmask = BPMDEF_NEAREDGE
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_NEAREDGE

                # Copy TAPEBUMP BPM bits to image mask
                bitmask = BPMDEF_TAPEBUMP
                mark = (bpm_im.mask & bitmask) != 0
                image.mask[mark] |= BADPIX_TAPEBUMP

                # Mark correctable pixels.
                # Pixels flagged as BPMDEF_BIAS_COL and BPMDEF_FUNKY_COL may be correctable.
                # We flag them in the image as bad (BADPIX_BPM), but - if fix_columns is run,
                # the BADPIX_BPM flag will be cleared and the BADPIX_FIX flag will be set
                # For each column find the number of pixels flagged as BIAS_HOT and BIAS_COL
                N_BIAS_HOT = np.sum((bpm_im.mask & BPMDEF_BIAS_HOT) > 0, axis=0)
                N_BIAS_COL = np.sum((bpm_im.mask & BPMDEF_BIAS_COL) > 0, axis=0)
                maskwidth=bpm_im.mask.shape[1]
                # First do columns with N_BIAS_COL set for 1 or more pixels
                biascols=np.arange(maskwidth)[(N_BIAS_COL > 0)]
                for icol in biascols:
                  #Clear FUNKY_COL bit if set for all pixels in this column
                  #The reason for clearing the bit is that the FUNKY_COL detection is
                  #sensitive to hot bias pixels and may flag those columns by "mistake"
                  #First clear BAD BPM bit if set because of funky column
                  image.mask[:,icol][bpm_im.mask[:,icol]==BPMDEF_FUNKY_COL] &= ~BADPIX_BPM
                  bpm_im.mask[:,icol] -= (bpm_im.mask[:,icol] & BPMDEF_FUNKY_COL )
                  #Correctable columns have exactly 1 BIAS_HOT pixel
                  if N_BIAS_HOT[icol] == 1:
                    #Correctable pixels have BIAS_COL bit set
                    bpm_im.mask[:,icol][(bpm_im.mask[:,icol]&BPMDEF_BIAS_COL)>0] |= BPMDEF_CORR
                    logger.info('Column '+str(icol)+' has 1 hot pixel and is correctable.')
                  else:
                    logger.info('Column '+str(icol)+' has '+str(N_BIAS_HOT[icol])+' hot pixels and is NOT correctable.')

                #Now do columns with FUNKY_COL set.  Note that the FUNKY_COL bits have been cleared above
                #for hot bias columns
                N_FUNKY_COL = np.sum((bpm_im.mask & BPMDEF_FUNKY_COL) > 0, axis=0)
                funkycols=np.arange(maskwidth)[(N_FUNKY_COL > 0)]
                for icol in funkycols:
                  #Correctable pixels have FUNKY_COL bit set
                  bpm_im.mask[:,icol][(bpm_im.mask[:,icol]&BPMDEF_FUNKY_COL)>0] |= BPMDEF_CORR
                  logger.info('Column '+str(icol)+' is funky and correctable.')

 
                image[kw] = time.asctime(time.localtime())
                image.write_key(kw, time.asctime(time.localtime()),
                                comment = 'Construct mask from BPM')
                if bpm_im.sourcefile is None:
                    image.write_key('BPMFIL', 'UNKNOWN', comment='BPM file used to build mask')
                else:
                    image.write_key('BPMFIL', path.basename(bpm_im.sourcefile), comment='BPM file used to build mask')
                        
                logger.debug('Finished applying BPM')

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
