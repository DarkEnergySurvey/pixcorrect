#!/usr/bin/env python
"""Produce a miniature sky-level image of a single CCD by taking
medians of boxes in image
"""

import ctypes
from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESDataImage, DESImage
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo

# Which section of the config file to read for this step
config_section = 'skycompress'

class SkyCompress(PixCorrectImStep):
    description = "Produce compressed image of sky background"
    step_name = config_section
    propagate = ['FILTER','DATE-OBS','EXPNUM','CCDNUM','DETPOS','INSTRUME',
                 'BAND','NITE']
    # Header keywords to copy from source image into compressed image
    
    @classmethod
    def __call__(cls, image, skyfilename, blocksize, bitmask):
        """Produce compressed image of sky background

        :Parameters:
            - `image`: the DESImage to be compressed.  
            - `skyfilename`: filename for the output compressed sky image
            - `blocksize`: side length of squares in which medians are taken
            - `bitmask`: Bitmask that will be or'ed with mask plane of image (if
                          any) to mark pixels to be ignored in calculating block
                          median.
        """
 
        logger.info('Compressing sky')

        nx = image.data.shape[1] / blocksize
        ny = image.data.shape[0] / blocksize
        ## ??? Check that blocksize evenly divides data.shape, else the following
        ## reshape will fail

        # Apply bit mask to the mask plane if any.  Superpixels
        # with no unmasked pixels will be filled with value -1
        if image.mask is None:
            sky = np.median(image.data.reshape(ny,blocksize,nx,blocksize)\
                            .swapaxes(1,2)\
                            .reshape(ny,nx,blocksize*blocksize), axis=2)
        else:
            data = np.ma.array(image.data, mask= (image.mask & bitmask),
                                fill_value=-1.)
            sky = np.ma.median(data.reshape(ny,blocksize,nx,blocksize)\
                               .swapaxes(1,2)\
                               .reshape(ny,nx,blocksize*blocksize), axis=2)
            sky = np.ma.getdata(sky)
        
        # Create HDU for output image, add some header info, save output to file
        outimage = DESDataImage(sky)
        outimage['BLOCKSIZ'] = blocksize
        outimage.copy_header_info(image, cls.propagate, require=False)
        ## ?? catch exception from write error below?
        outimage.save(skyfilename)

        logger.debug('Finished sky compression')
        ret_code=0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for sky compression:

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        ### ?? Put defaults in here ??

        if config.has_option(cls.step_name,'blocksize'):
            blocksize = config.getint(cls.step_name, 'blocksize')
        else:
            blocksize = skyinfo.DEFAULT_BLOCKSIZE
        if config.has_option(cls.step_name,'bitmask'):
            bitmask = config.getint(cls.step_name, 'bitmask')
        else:
            bitmask = skyinfo.DEFAULT_SKYMASK
        skyfilename = config.get(cls.step_name, 'skyfilename')
        logger.info('Sky compression will be done for %s' % image)
    
        ret_code = cls.__call__(image, skyfilename, blocksize, bitmask)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--skyfilename', default='minisky.fits',
                            help='Filename for compressed sky image')
        parser.add_argument('--blocksize', type=int, default=skyinfo.DEFAULT_BLOCKSIZE,
                            help='Size of squares in which median is taken for sky')
        parser.add_argument('--bitmask', type=int, default=skyinfo.DEFAULT_SKYMASK,
                            help='Mask image bits for pixels to ignore in sky estimate')
        return

    @classmethod
    def run(cls, config):
        """Execute the sky compress step, opening the input FITS file first.
        Need to override the PixCorrectImStep version since this step has no
        output version of the image.
        """
        in_fname = config.get(cls.step_name, 'in')
        try:
            ccdnum = config.getint(cls.step_name, 'ccdnum')
            image = DESImage.load(in_fname, ccdnum=ccdnum)
        except NoOptionError:
            image = DESImage.load(in_fname)

        ret_code = cls.step_run(image, config)
        return ret_code

sky_compress = SkyCompress()

# internal functions & classes

if __name__ == '__main__':
    sky_compress.main()
