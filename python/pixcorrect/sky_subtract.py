#!/usr/bin/env python
"""
Fit templates to a mini-sky image of a full exposure, returns best-fit coefficients
and statistics on the residuals to fit.
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
    description = "Subtract sky from images based on principal-component fit"
    step_name = config_section
    
    @classmethod
    def __call__(cls, image, fit_filename, pc_filename,
                sky_weight, domeflat, bitmask):
        """
        Subtract sky from image using previous principal-components fit. Optionally
        build weight image from fitted sky, in which case the dome flat is needed and
        proper gain values are expected in the image header.

        :Parameters:
            - `image`: DESImage that has been flattened with dome already and fit
            - `fit_filename`: filename with the coefficients from minisky fitting
            - `pc_filename`: filename for the stored full-res sky principal components
            - `sky_weight`: set True to (re-)build the weight image.
            - `domeflat`: object whose data property is the 2d dome flat image, needed
                          for building the weight image.
            - `bitmask`: which bits in the mask image signal that weight=0.
        """
 
        logger.info('Subtracting sky')

        mini = skyinfo.MiniDecam.load(fit_filename)
        templates = skyinfo.SkyPC.load(pc_filename)
        sky = templates.sky(mini.coeffs)
        image.data -= sky
        logger.debug('Finished sky subtraction')
        # ?? Any special treatment of masked pixels?

        if sky_weight:
            logging.info('Constructing weight image from sky image')
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
            """

            # Check for needed input
            if domeflat is None:
                raise SkyError('sky_fit needs dome flat when sky_weight=True')

            # Transform the sky image into a variance image
            for amp in decaminfo.amps:
                sec = scan_fits_section(image,'DATASEC'+amp)
                sec = (slice(sec[0]-1,sec[1]),slice(sec[2]-1,sec[3]))
                sky[sec] *= image['DOMEMED'+amp]/image['GAIN'+amp]
                sky[sec] += (image['RDNOISE'+amp]*image['DOMEMED'+amp])**2 / dome[sec]
            
            sky /= dome
            # Weight is to be inverse variance:
            image.weight = 1./sky
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
            clip_sigma = config.getint(cls.step_name, 'bitmask')
        else:
            clip_sigma = skyinfo.DEFAULT_BITMASK

        fit_filename = config.get(cls.step_name, 'fitfilename')
        pc_filename = config.get(cls.step_name, 'pcfilename')
        if config.has_option(cls.step_name,'skyweight'):
            sky_weight = config.getboolean(cls.step_name, 'skyweight')
        else:
            sky_weight = skyinfo.DEFAULT_SKY_WEIGHT
        if sky_weight:
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
        parser.add_argument('--skyweight', type=boolean, default=skyinfo.DEFAULT_SKY_WEIGHT,
                            help='Construct weight image from fitted sky')
        parser.add_argument('--domefilename',type=str,
                            help='Filename for dome flat (for skyweight=True)')
        parser.add_argument('--bitmask',type=int,default=skyinfo.DEFAULT_BITMASK,
                            help='which bits in mask image imply weight=0')
        return

sky_subtract = SkySubtract()

# internal functions & classes

if __name__ == '__main__':
    sky_subtract.main()
