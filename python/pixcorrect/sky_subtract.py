#!/usr/bin/env python
"""
Subtract sky from image by summing sky principal components with pre-computed coefficients
for this exposure.
Also create or overwrite a weight image, using read noise, flat-field noise, and shot noise either from all counts or just from sky.
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESDataImage, DESImage, weight_dtype, section2slice
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo
from pixcorrect.skyinfo import SkyError
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'skysubtract'

class SkySubtract(PixCorrectImStep):
    description = "Subtract sky from images based on principal-component fit and calculate" +\
      " weight image"
    
    step_name = config_section
    
    @classmethod
    def __call__(cls, image, fit_filename, pc_filename,
                weight, dome):
        """
        Subtract sky from image using previous principal-components fit. Optionally
        build weight image from fitted sky or all counts, in which case the dome flat
        is needed and proper gain values are expected in the image header.

        :Parameters:
            - `image`: DESImage that has been flattened with dome already and fit
            - `fit_filename`: filename with the coefficients from minisky fitting.  Sky
                              subtraction is skipped if this is None.
            - `pc_filename`: filename for the stored full-res sky principal components
            - `weight`: 'none' to skip weights, 'sky' to calculate weight at sky level,
                         'all' to use all counts
            - `dome`: DESImage for the dome flat, needed if weight is not 'none'.
        """
 
        if weight=='sky' and fit_filename is None:
            raise SkyError('Cannot make sky-only weight map without doing sky subtraction')
        
        if fit_filename is not None:
            logger.info('Subtracting sky')
            mini = skyinfo.MiniDecam.load(fit_filename)
            templates = skyinfo.SkyPC.load(pc_filename)
            sky = templates.sky(mini.coeffs)
            image.data -= sky
            logger.debug('Finished sky subtraction')

        if weight=='none':
            do_weight = False
            sky_weight = False
        elif weight=='sky':
            do_weight = True
            sky_weight = True
        elif weight=='all':
            do_weight = True
            sky_weight = False
        else:
            raise SkyError('Invalid weight value: ' + weight)

        if do_weight:
            if dome is None:
                raise SkyError('sky_subtract needs dome flat when making weights')
        
            if sky_weight:
                logger.info('Constructing weight image from sky image')
                data = sky
            else:
                logger.info('Constructing weight image from all counts')
                data = image.data + sky

            if image.weight is not None or image.variance is not None:
                image.weight = None
                image.variance = None
                logger.warning('Overwriting existing weight image')
                
            """
            We assume in constructing the weight (=inverse variance) image that
            the input image here has been divided by the dome flat already, and that
            its GAIN[AB] keywords are correct for a pixel that has been divided
            by the FLATMED[AB] of the flat image.  So the number of *electrons* that
            were read in a pixel whose current value=sky is
            e = sky * (dome/FLATMED) * GAIN


            The variance has three parts: read noise, and sky Poisson noise, and
            multiplicative errors from noise in the flat field.
            The read noise variance, in electrons, is
            Var = RDNOISE^2
            ...and the shot noise from sky was, in electrons,
            Var = sky * (dome/FLATMED) * GAIN

            This means the total variance in the image, in its present form, is

            Var = (RDNOISE * FLATMED / dome / GAIN)^2 + (FLATMED/GAIN)*sky/dome

            We can also add the uncertainty propagated from shot noise in the dome flat,
            if the dome image has a weight or variance.  In which case we would add

            Var += var(dome) * sky^2 / dome^2

            (remembering that sky has already been divided by the dome).

            If sky_weight = False, we can substitute the image data for sky in the above
            calculations.
            """

            # Transform the sky image into a variance image
            var = np.array(data, dtype = weight_dtype)
            for amp in decaminfo.amps:
                sec = section2slice(image['DATASEC'+amp])
                invgain = (image['FLATMED'+amp]/image['GAIN'+amp]) / dome.data[sec]
                var[sec] += image['RDNOISE'+amp]**2 * invgain
                var[sec] *= invgain
            # Add noise from the dome flat shot noise, if present
            if dome.weight is not None:
                var += data * data / (dome.weight*dome.data * dome.data)
            elif dome.variance is not None:
                var += data * data * dome.variance / (dome.data * dome.data)
                
            image.variance = var

            ## ??? Add SKYVAR[AB] keywords for Mangle ???
            
            logger.debug('Finished weight construction')

        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for sky subtraction

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        if config.has_option(cls.step_name,'fitfilename'):
            fit_filename = config.get(cls.step_name, 'fitfilename')
        else:
            fit_filename = None

        if config.has_option(cls.step_name,'pcfilename'):
            pc_filename = config.get(cls.step_name, 'pcfilename')
        else:
            pc_filename = None

        weight = config.get(cls.step_name,'weight')

        if config.has_option(cls.step_name,'domefilename'):
            dome_filename = config.get(cls.step_name,'domefilename')
            dome = DESImage.load(dome_filename)
        else:
            dome = None

        logger.info('Sky fitting output to %s' % image)
    
        ret_code = cls.__call__(image, fit_filename, pc_filename,
                                weight, dome)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--fitfilename',type=str,
                            help='Filename for minisky FITS image with PC coefficients')
        parser.add_argument('--pcfilename',type=str,
                            help='Filename for full-res sky principal components')
        parser.add_argument('--domefilename',type=str,
                            help='Filename for dome flat (for weight calculation)')
        parser.add_argument('--weight', choices=('sky','all','none'),
                            default=skyinfo.DEFAULT_WEIGHT,
                            help='Construct weight from sky photons, ' \
                                 'from all photons, or not at all')
        return

sky_subtract = SkySubtract()

# internal functions & classes

if __name__ == '__main__':
    sky_subtract.main()
