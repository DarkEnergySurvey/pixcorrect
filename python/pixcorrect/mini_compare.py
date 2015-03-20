#!/usr/bin/env python
"""
Compare two compressed DES images, report size of deviations
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESDataImage, DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo

# Which section of the config file to read for this step
config_section = 'minicompare'

class MiniCompare(PixCorrectImStep):
    description = "Compare input compressed DECam image to reference after matching normalizations." \
      "The header of the output image will also contain keywords RMS and WORST giving the rms" \
      " and maximum fractional deviations, and a keyword FACTOR giving the overall flux factor" \
      " used to normalize them."
    step_name = config_section
    
    @classmethod
    def __call__(cls, in_filename, ref_filename, out_filename):
        """
        Compare two compressed DES images, report size of deviations.

        :Parameters:
            - `in_filename`: the compressed DES image ("mini-sky" format) to check
            - `ref_filename`: the reference image to compare to, must have same compression as in_filename
            - `out_filename`: output image showing fractional residuals of input to reference after matching
                              normalizations.  The header of the output image will also contain keywords RMS
                              and WORST giving the rms and maximum fractional deviations, and a keyword FACTOR
                              giving the overall flux factor used to normalize them.
        """
 
        logger.info('Comparing compressed image' + in_filename + " to reference " + ref_filename)

        indata = skyinfo.MiniDecam.load(in_filename)
        ref = skyinfo.MiniDecam.load(ref_filename)

        if indata.blocksize != ref.blocksize or indata.invalid != ref.invalid or indata.halfS7 != ref.halfS7:
            raise SkyError("Input and reference are not matching compressions of DECam")

        resid = indata.vector() / ref.vector()
        factor = np.median(resid)
        resid /= factor
        resid -= 1.
        rms = np.std(resid)
        worst = np.max(np.abs(resid))
        indata.fill_from(resid)
        indata.header['FACTOR'] = factor
        indata.header['RMS'] = rms
        indata.header['WORST'] = worst
        indata.save(out_filename)

        logger.info('Normalization factor: %f' % factor)
        logger.info('RMS deviation: %f' % rms)
        logger.info('Worst deviation: %f' % worst)

        # Create a one-line binary fits table to hold the coefficients
        logger.debug('Finished image comparison')
        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, config):
        """Customized execution for img comparison. No input image nor output

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        in_filename = config.get(cls.step_name, 'in')
        out_filename = config.get(cls.step_name, 'out')
        ref_filename = config.get(cls.step_name, 'ref')

        ret_code = cls.__call__(in_filename, ref_filename, out_filename)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('-r','--ref',type=str,
                            help='Filename for reference compressed image')
        return

    @classmethod
    def run(cls, config):
        """Execute the image comparison step.
        Need to override the PixCorrectImStep version since this step does not
        have full DECam images as input nor output.
        """
        ret_code = cls.step_run(config)
        return ret_code

mini_compare = MiniCompare()

# internal functions & classes

if __name__ == '__main__':
    mini_compare.main()
