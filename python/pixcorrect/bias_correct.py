#!/usr/bin/env python
"""Apply a bias correction to a raw DES image 
"""

from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once, items_must_match
from despyfits.DESImage import DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'bias'


class BiasCorrect(PixCorrectImStep):
    description = "Apply a bias correction to an image"
    step_name = config_section

    @classmethod
    @do_once(1,'DESBIAS')
    def __call__(cls, image, bias_im):
        """Apply a bias correction to an image

        :Parameters:
            - `image`: the DESImage to apply a bias correction
            - `bias_im`:  the bias correction image to apply

        Applies the correction "in place." Also creates BAND and NITE
        keywords if they are not present.
        """
 
        logger.info('Applying Bias')
        # Check that bias and data are from same CCD
        try:
            items_must_match(image, bias_im, 'CCDNUM')
        except:
            return 1
        image.data -= bias_im.data
        # If we have two weight images, add variance of the bias to the image's
        if (image.weight is not None or image.variance is not None):
            if bias_im.weight is not None:
                var = image.get_variance()
                var += 1./bias_im.weight
            elif bias_im.variance is not None:
                var = image.get_variance()
                var += bias_im.variance
        logger.debug('Finished applying Bias')
        if bias_im.sourcefile is None:
            image.write_key('BIASFIL', 'UNKNOWN', comment='Bias correction file')
        else:
            image.write_key('BIASFIL', path.basename(bias_im.sourcefile), comment='Bias correction file')
        # Also create the BAND and NITE keywords if they are not present
        try:
            image['BAND']
        except:
            image['BAND'] = decaminfo.get_band(image['FILTER'])
        try:
            image['NITE']
        except:
            image['NITE'] = decaminfo.get_nite(image['DATE-OBS'])
            
        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `bias`: the bias image to apply

        """

        bias_fname = config.get(cls.step_name, 'bias')
        logger.info('reading Bias from %s'% bias_fname)
        bias_im = DESImage.load(bias_fname)
    
        ret_code = cls.__call__(image, bias_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the bias correction
        """
        parser.add_argument('--bias', nargs=1, default=None,
                            help='Bias correction image')

bias_correct = BiasCorrect()

# internal functions & classes

if __name__ == '__main__':
    bias_correct.main()
