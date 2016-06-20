#!/usr/bin/env python
"""
Fit templates to a mini-sky image of a full exposure, returns best-fit coefficients
and statistics on the residuals to fit.
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger, items_must_match
from despyfits.DESImage import DESDataImage, DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo

# Which section of the config file to read for this step
config_section = 'skyfit'

class SkyFit(PixCorrectImStep):
    description = "Fit coefficients of sky templates to mini-sky image"
    step_name = config_section
    
    @classmethod
    def __call__(cls, in_filename, out_filename, pc_filename, clip_sigma):
        """
        Fit coefficients of sky templates to mini-sky FITS image from an exposure.

        :Parameters:
            - `in_filename`: filename for exposure mini-sky image
            - `out_filename`: filename for the output coefficients and residual sky data
            - `pc_filename`: filename for the stored mini-sky principal component array
            - `clip_sigma`: Number of sigma to mark outliers ignored from fitting & stats
        """
 
        logger.info('Fitting sky')

        mini = skyinfo.MiniDecam.load(in_filename)
        templates = skyinfo.MiniskyPC.load(pc_filename)
        try:
            # Insure using the correct filter's PCA
            items_must_match(mini.header,templates.header,'BAND')
        except:
            return 1
        templates.fit(mini, clip_sigma)
        mini.save(out_filename)

        # Create a one-line binary fits table to hold the coefficients
        logger.debug('Finished sky fitting')
        ret_code=0
        return ret_code

    @classmethod
    def step_run(cls, config):
        """Customized execution for sky combination.  Note there is NO input image nor output

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        if config.has_option(cls.step_name,'clipsigma'):
            clip_sigma = config.getfloat(cls.step_name, 'clipsigma')
        else:
            clip_sigma = skyinfo.DEFAULT_CLIP_SIGMA

        in_filename = config.get(cls.step_name, 'infilename')
        out_filename = config.get(cls.step_name, 'outfilename')
        pc_filename = config.get(cls.step_name, 'pcfilename')

        logger.info('Sky fitting output to %s' % out_filename)
    
        ret_code = cls.__call__(in_filename, out_filename, pc_filename, clip_sigma)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--infilename',type=str,
                            help='Filename for input minisky FITS image to fit')
        parser.add_argument('--outfilename',type=str,
                            help='Filename for minisky FITS image with fit results/resids')
        parser.add_argument('--pcfilename',type=str,
                            help='Filename for minisky principal components')
        parser.add_argument('--clip_sigma', type=float, default=skyinfo.DEFAULT_CLIP_SIGMA,
                            help='Rejection threshold for robust fitting/statistics')
        return

    @classmethod
    def run(cls, config):
        """Execute the sky combine step.
        Need to override the PixCorrectImStep version since this step has no
        input nor output image
        """
        ret_code = cls.step_run(config)
        return ret_code

sky_fit = SkyFit()

# internal functions & classes

if __name__ == '__main__':
    sky_fit.main()
