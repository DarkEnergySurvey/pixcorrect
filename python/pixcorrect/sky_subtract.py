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
from despyfits.DESImage import DESDataImage, DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo

# Which section of the config file to read for this step
config_section = 'skysubtract'

class SkySubtract(PixCorrectImStep):
    description = "Subtract sky from images based on principal-component fit and calculate" +\
      " weight image"
    
    step_name = config_section
    
    @classmethod
    def __call__(cls, image, fit_filename, pc_filename,
                weight, domeflat, bitmask):
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
            - `domeflat`: DESImage for the dome flat, needed if weight=True.
            - `bitmask`: which bits in the mask image signal that weight=0.
        """
 
        if weight=='sky' and fit_filename is None:
            raise SkyError('Cannot make sky-only weight map without doing sky subtraction')
        
        if domeflat is None:
            raise SkyError('sky_fit needs dome flat when weight=True')
        
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
        elif weight=='sky':
            do_weight = True
            sky_weight = False
        else:
            raise SkyError('Invalid weight value: ' + weight)

        if do_weight:
            if sky_weight:
                logging.info('Constructing weight image from sky image')
                data = sky
            else:
                logging.info('Constructing weight image from all counts')
                data = image.data

            # ?? change to: if image.weight is not None or image.variance is not None:
            if image.weight is not None:
                logging.warning('Overwriting existing weight image')
                
            """
            We assume in constructing the weight (=inverse variance) image that
            the input image here has been divided by domeflat already, and that
            its RDNOISE[AB] and GAIN[AB] are correct for a pixel that has been divided
            by the DOMEMED[AB] of the dome flat image.

            Also the RDNOISE and GAIN reflect whatever re-scaling has been done - so
            if the image was *divided* by SCALE, we changed
            RDNOISE -> RDNOISE / SCALE
            GAIN -> GAIN * SCALE
            
            The variance has two parts: read noise and sky Poisson.
            The read noise variance, before scaling by dome, was
            Var = RN^2 = (RDNOISE*DOMEMED)^2
            ...and the shot noise from sky was
            Var = (original sky) / (original gain)
                = (sky * dome) / (GAIN / DOMEMED)
            This means the total variance in the image, after dividing by dome, is

            Var = (RDNOISE * DOMEMED / dome)^2 + (DOMEMED/GAIN)*sky/dome

            We can also add the uncertainty propagated from shot noise in the dome flat,
            if the dome image has a weight or variance.  In which case we would add

            Var += var(dome) * sky^2 / dome^2

            (remembering that sky has already been divided by the dome).

            If sky_weight = False, we can substitute the image data for sky in the above
            calculations.
            """

            # Transform the sky image into a variance image
            var = np.array(data, dtype = DESImage.weight_dtype)
            for amp in decaminfo.amps:
                sec = section2slice(image['DATASEC'+amp])
                var[sec] *= image['FLATMED'+amp]/image['GAIN'+amp]
                var[sec] += (image['RDNOISE'+amp]*image['FLATMED'+amp])**2 / dome.data[sec]
            var /= dome.data
            # Add noise from the dome flat shot noise, if present
            # ??? change when variance is available
            if dome.weight is not None:
                var += data * data * / (dome.weight*dome.data * dome.data)
                
            # ?? Change when variance is available:
            # Weight is to be inverse variance:
            image.weight = 1./var
            # Apply a mask to the weight image
            image.weight[ image.mask | bitmask] = 0.

            logger.debug('Finished weight construction')

        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, image, config):
        """Customized execution for sky subtraction

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        if config.has_option(cls.step_name,'bitmask'):
            bitmask = config.getint(cls.step_name, 'bitmask')
        else:
            bitmask = skyinfo.DEFAULT_BITMASK

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
            domeflat = DESDataImage.load(dome_filename)
        else:
            domeflat = None

        logger.info('Sky fitting output to %s' % image)
    
        ret_code = cls.__call__(image, fit_filename, pc_filename,
                                sky_weight, domeflat, bitmask)
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
        parser.add_argument('--bitmask',type=int,default=skyinfo.DEFAULT_BITMASK,
                            help='which bits in mask image imply weight=0')
        parser.add_argument('--weight', choices=('sky','all','none'),
                            default=skyinfo.DEFAULT_WEIGHT,
                            help='Construct weight from sky photons, ' \
                                 'from all photons, or not at all')
        return

sky_subtract = SkySubtract()

# internal functions & classes

if __name__ == '__main__':
    sky_subtract.main()
