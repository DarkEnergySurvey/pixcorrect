#!/usr/bin/env python
"""Apply a bias correction to a raw DES image 
"""

from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once
from despyfits.DESImage import DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'imgdiff'


class ImageDiff(PixCorrectImStep):
    description = "Take difference between input image and another `comparison` image"
    step_name = config_section

    @classmethod
    def __call__(cls, image, comp_im):
        """Take difference between image and another `comparison` image "

        :Parameters:
            - `image`: the DESImage to apply a bias correction
            - `comp_im`:  comparison image to be subtracted

        Applies the correction "in place"
        """
 
        logger.info('Taking Difference')
        image.data -= comp_im.data
        # If we have two weight images, add variance of the bias to the image's
        if (image.weight is not None or image.variance is not None):
            if comp_im.weight is not None:
                var = image.get_variance()
                var += 1./comp_im.weight
            elif comp_im.variance is not None:
                var = image.get_variance()
                var += comp_im.variance
        logger.debug('Finished taking difference')
        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for taking difference between an image and a comparison

        :Parameters:
            - `image`: the DESImage on which to operate
            - `comp`: the comparison image (to be subtracted)

        """

        comp_fname = config.get(cls.step_name, 'comp')
        logger.info('reading Comparison image from %s'% comp_fname)
        comp_im = DESImage.load(comp_fname)
    
        ret_code = cls.__call__(image, comp_im)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to comparison subtraction
        """
        parser.add_argument('--comp', nargs=1, default=None,
                            help='Comparison image (to be subtracted)')

image_diff = ImageDiff()

# internal functions & classes

if __name__ == '__main__':
    image_diff.main()
