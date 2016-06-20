#!/usr/bin/env python
"""Apply a flat correction to a raw DES image 
"""

from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, do_once, items_must_match
#from despyfits.DESImage import DESImage
from despyfits.DESImage import DESImage, section2slice 
from despyfits import maskbits
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'flat'

class FlatCorrect(PixCorrectImStep):
    description = "Apply a flat field correction to an image"
    step_name = config_section

    @classmethod
    def _doit(cls, image, flat_im):
        """Apply a flat field correction to an image - used for both dome
        and star flats.

        :Parameters:
            - `image`: the DESImage to apply a bias correction
            - `flat_im`:  the flat correction image to apply

        Applies the correction "in place"
        """
        logger.info('Applying Flat')
        
        # Check that flat and data are from same CCD and filter
        try:
            image['BAND']
        except:
            # Give image a BAND from its FILTER if it's not there
            image['BAND'] = decaminfo.get_band(image['FILTER'])
        try:
            items_must_match(image, flat_im, 'CCDNUM','BAND')
        except:
            return 1
        
        # Apply flat to the data
        image.data /= flat_im.data

        # Update variance or weight image if it exists
        if image.weight is not None:
            image.weight *= flat_im.data*flat_im.data
        if image.variance is not None:
            image.variance /= flat_im.data*flat_im.data

        # If mask image exists, mark as BADPIX_BPM any pixels that have
        # non-positive flat and are not already flagged.
        if image.mask is not None:
            # Find flat-field pixels that are invalid but not already bad for
            # one of these reasons:
            badmask = maskbits.BADPIX_BPM +\
              maskbits.BADPIX_BADAMP +\
              maskbits.BADPIX_EDGE
            badflat = np.logical_and( flat_im.data <= 0.,
                                      image.mask & badmask)
            mark_these = np.where(badflat.flatten())[0]
            image.mask.flatten()[mark_these] |= maskbits.BADPIX_BPM
            
        # If a weight or variance image already exists, add to it any additional
        # variance from the flat:
        if (image.weight is not None or image.variance is not None):
            if flat_im.weight is not None:
                var = image.get_variance()
                f2 = flat_im.data * flat_im.data
                var *= f2
                var += image.data*image.data/(flat_im.weight*f2)
            elif flat_im.variance is not None:
                var = image.get_variance()
                f2 = flat_im.data * flat_im.data
                var *= f2
                var += image.data*image.data*flat_im.variance/f2

        # Update header keywords for rescaling
        saturate = 0.
        scales = []
        for amp in decaminfo.amps:
            # Acquire the typical scaling factor for each amp from the flat
            scalekw = 'FLATMED'+amp
            if scalekw in flat_im.header.keys():
                # Already stored in the flat's header:
                scale = flat_im[scalekw]
            else:
                # Figure it out ourselves from median of a subsample:
#                sec = DESImage.section2slice(image['DATASEC'+amp])
                sec = section2slice(image['DATASEC'+amp])
                scale = np.median(flat_im.data[sec][::4,::4])
            scales.append(scale)
            if scalekw in image.header.keys():
                # Add current scaling to any previous ones
                image[scalekw] = image[scalekw]*scale
            else:
                image[scalekw] = scale
            image['GAIN'+amp] = image['GAIN'+amp] * scale
            image['SATURAT'+amp] = image['SATURAT'+amp] / scale
            # Scale the SKYVAR if it's already here
            kw = 'SKYVAR'+amp
            if kw in image.header.keys():
                image[kw] = image[kw] / (scale*scale)
            saturate = max(saturate, image['SATURAT'+amp])
        # The SATURATE keyword is assigned to maximum of the amps' values.
        image['SATURATE'] = saturate
            
        # Some other keywords that we will adjust crudely with mean rescaling
        # if they are present:
        scale = np.mean(scales)
        for kw in ('SKYBRITE','SKYSIGMA'):
            if kw in image.header.keys():
                image[kw] = image[kw] / scale

        logger.debug('Finished applying Flat')
        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `flat`: the bias image to apply

        """

        flat_fname = config.get(cls.step_name, 'flat')
        logger.info('Reading flat correction from %s'% flat_fname)
        flat_im = DESImage.load(flat_fname)
    
        ret_code = cls.__call__(image, flat_im)
        return ret_code

    @classmethod
    @do_once(1,'DESFLAT')
    def __call__(cls, image, flat_im):
        """Apply a flat field correction to an image

        :Parameters:
            - `image`: the DESImage to flatten
            - `flat_im`:  the flat correction image to apply

        Applies the correction "in place"
        """
        ret_code = cls._doit(image, flat_im)

        if flat_im.sourcefile is None:
            image.write_key('FLATFIL', 'UNKNOWN', comment='Dome flat correction file')
        else:
            image.write_key('FLATFIL', path.basename(flat_im.sourcefile), comment='Dome flat correction file')

        return ret_code
            
    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--flat', nargs=1, default=None,
                            help='Flat field correction image')

flat_correct = FlatCorrect()

# internal functions & classes

if __name__ == '__main__':
    flat_correct.main()
