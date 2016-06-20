#!/usr/bin/env python
"""Apply brighter-fatter correction to an image
"""

from os import path
import numpy as np
from pixcorrect.corr_util import logger, do_once
from despyfits.DESImage import DESImage, section2slice, data_dtype
from despyfits.maskbits import *
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo
from pixcorrect.bfinfo import *

# Which section of the config file to read for this step
config_section = 'bf'

class BFCorrect(PixCorrectImStep):
    description = "Apply brighter-fatter correction to an image"
    step_name = config_section

    @classmethod
    @do_once(1,'DESBFC')
    def __call__(cls, image, bffile, bfmask):
        """
        Apply brighter-fatter correction to an image, and set BADPIX_SUSPECT bit
        in mask image for pixels adjacent to those with unknown collected charge.

        :Parameters:
            - `image`: the DESImage to operate upon.  Must have mask plane present.
            - `bffile`: name of FITS file holding the brighter-fatter coefficients
            - `bfmask`: which bits in the mask will trigger pixel being ignored - should
                        denote those pixels having unknown amount of charge during integration.
        """

        logger.info('Start brighter-fatter correction')

        if image.mask is None:
            raise BFError("Missing mask image for bf_correct")

        detpos = image['DETPOS'].strip()
        logger.info('reading BF corrections from %s' % bffile)
        bf = BFKernel(bffile, detpos)

        ignore = np.logical_or(np.isinf(image.data), np.isnan(image.data))
        ignore = np.logical_or(ignore, (image.mask & bfmask) != 0)
            
        # Get a median sky level and replace bad pixels with it when deriving kernel
        # Also put image into electron units, if not already.
        data = np.array(image.data)    
        for amp in decaminfo.amps:
            gain = image['GAIN'+amp]
            sec = section2slice(image['DATASEC'+amp])
            if gain != 1:
                data[sec] *= gain
            sky = np.median(data[sec][::4,::4])
            data[sec][ignore[sec]] = sky

        # Convolve data with R kernel to get right-hand pixel shifts
        df = np.fft.rfft2(data)
        kernel = bf.kernelR(data.shape)
        shift = np.fft.irfft2(df * np.fft.rfft2(kernel))
        # Multiply by border charge to get amount of charge to move.
        charge = 0.5*(data[:,:-1] + data[:,1:])*shift[:,:-1]
        # Do not shift charge into or out of bad pixels
        charge[ignore[:,:-1]] = 0.
        charge[ignore[:,1:]] = 0.

        # Adjust data for this shift
        out = np.array(image.data)
        for amp in decaminfo.amps:
            # Redo the temporary gain correction:
            gain = image['GAIN'+amp]
            sec = section2slice(image['DATASEC'+amp])
            if gain != 1:
                out[sec] *= gain
        out[:,1:] -= charge
        out[:,:-1] += charge

        # Now do the upper-edge pixel shifts  ??? Add gain factor here & T?
        kernel = bf.kernelT(data.shape)
        shift = np.fft.irfft2(df * np.fft.rfft2(kernel))
        # Multiply by border charge to get amount of charge to move.
        charge = 0.5*(data[:-1,:] + data[1:,:]) * shift[:-1,:]
        # Do not shift charge into or out of bad pixels
        charge[ignore[:-1,:]] = 0.
        charge[ignore[1:,:]] = 0.
        # Adjust data for this shift
        out[1:,:] -= charge
        out[:-1,:] += charge

        # Undo the gain correction if we made it originally:
        for amp in decaminfo.amps:
            gain = image['GAIN'+amp]
            sec = section2slice(image['DATASEC'+amp])
            if gain != 1.:
                out[sec] /= gain
        image.data = out

        # Set the SUSPECT flag for all pixels that were adjacent to
        # ignored pixels, as their b/f correction is off
        change_mask = np.zeros(image.mask.shape,dtype=bool)
        change_mask[:-1,:] |= ignore[1:,:]   # mask below
        change_mask[1:,:] |= ignore[:-1,:]   # mask above
        change_mask[:,:-1] |= ignore[:,1:]   # mask to left
        change_mask[:,1:] |= ignore[:,:-1]   # mask to right
        change_mask[ignore] = False          # Don't mask what's already bad
        image.mask[change_mask] |= BADPIX_SUSPECT

        image.write_key('BFCFIL', path.basename(bffile), comment='Brighter/fatter correction file')
        
        ret_code = 0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for B/F correction

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        if config.has_option(cls.step_name, 'bfmask'):
            bfmask = parse_badpix_mask(config.get(cls.step_name, 'bfmask'))
        else:
            bfmask = parse_badpix_mask(DEFAULT_BFMASK)

        bffile = config.get(cls.step_name,'bffile')
            
        ret_code = cls.__call__(image, bffile, bfmask)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('--bffile', 
                            help='B/F coefficients filename')
        parser.add_argument('--bfmask', default=DEFAULT_BFMASK,
                            help='Bitmask for pixels to ignore in B/F correction')

bf_correct = BFCorrect()

# internal functions & classes

if __name__ == '__main__':
    bf_correct.main()
