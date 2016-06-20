#!/usr/bin/env python
"""
Combine all of the individual CCD's mini-sky images into one FITS image for whole array
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from despyfits.DESImage import DESDataImage, DESImage
from pixcorrect.PixCorrectDriver import PixCorrectDriver, filelist_to_list
from pixcorrect import skyinfo

# Which section of the config file to read for this step
config_section = 'skycombine'

class SkyCombine(PixCorrectDriver):
    description = "Combine sky images of all CCDs in one exposure"
    step_name = config_section
    propagate = ['FILTER','DATE-OBS','EXPNUM','INSTRUME','BAND','NITE']
    # Header keywords to copy from a single-chip image into the output image
    
    @classmethod
    def __call__(cls, in_filenames, out_filename, mask_value, invalid):
        """
        Produce compressed image of sky background for entire exposure by
        combining minisky FITS images for single CCDs.  Missing input minisky
        will only generate a warning.

        :Parameters:
            - `in_filenames`: list of filenames of single-chip sky images
            - `out_filename`: filename for the output combined sky image
            - `mask_value`: value inserted into sky pixels with no data
            - `invalid`: list of detpos values for CCDs to be left out of sky image
        """
 
        logger.info('Combining sky')

        out = None
        # Insert each input image into the output image
        for f in in_filenames:
            try:
                small = DESDataImage.load(f)
            except (ValueError,IOError) as v:
                # A missing minisky file is not fatal:
                logger.warning('SkyCombine could not load minisky '+f)
                continue
            if small['DETPOS'].strip() in invalid:
                # Skip any unwanted CCDs
                continue
            blocksize = small['BLOCKSIZ']
            if out is None:
                out = skyinfo.MiniDecam(blocksize, mask_value, invalid)
                out.copy_header_info(small, cls.propagate, require=False)
                
            if blocksize != out.blocksize:
                raise skyinfo.SkyError('Mismatched blocksizes for SkyCombine')
            out.fill(small.data,small['DETPOS'].strip())
            # Issue warnings for mismatches of data, but skip if
            # quantities are not in headers
            for k in ('BAND','NITE','EXPNUM'):
                try:
                    v1 = out[k]
                    v2 = small[k]
                    if not v1==v2:
                        logger.warning('Mismatched {:s} in file {:s}: {} vs {}'.\
                                       format(k,f,v1,v2))
                except:
                    pass

        out.save(out_filename)

        logger.debug('Finished sky combination')
        ret_code=0
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to sky compression
        """
        parser.add_argument('--ccdnums',type=str,default=skyinfo.DEFAULT_CCDNUMS,
                            help='Range(s) of ccdnums to combine')
        parser.add_argument('--miniskyfiles',type=str,default=skyinfo.DEFAULT_MINISKY_FILES,
                            help='Filename template for single-chip minisky images')
        parser.add_argument('--miniskylist',type=str,default=None, help='File containing a list of single-chip minisky images')
        parser.add_argument('-o','--outfilename',type=str,
                            help='Filename for combined minisky FITS image')
        parser.add_argument('--mask_value', type=float, default=skyinfo.DEFAULT_MASK_VALUE,
                            help='Value of pixels without valid sky information')
        parser.add_argument('--invalid', type=str, default=skyinfo.DEFAULT_IGNORE,
                            help='Value(s) of DETPOS to ignore in sky image')
        return

    @classmethod
    def run(cls, config):
        """Customized execution for sky combination.  Note there is NO input image nor output

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        if config.has_option(cls.step_name,'maskvalue'):
            mask_value = config.getfloat(cls.step_name, 'maskvalue')
        else:
            mask_value = skyinfo.DEFAULT_MASK_VALUE

        if config.has_option(cls.step_name,'invalid'):
            baddet = config.get(cls.step_name, 'invalid')
        else:
            baddet = skyinfo.DEFAULT_IGNORE
        invalid = baddet.split(',')

        if config.has_option(cls.step_name,'ccdnums'):
            ccdranges = config.get(cls.step_name,'ccdnums')
        else:
            ccdranges = skyinfo.DEFAULT_CCDNUMS
        ccdnumlist = skyinfo.parse_ranges(ccdranges)

        if config.has_option(cls.step_name, 'miniskylist'):
            miniskylist = config.get(cls.step_name,'miniskylist')
            in_filenames=filelist_to_list(miniskylist)
        else:
            if config.has_option(cls.step_name, 'miniskyfiles'):
                miniskyfiles = config.get(cls.step_name,'miniskyfiles')
            else:
                miniskyfiles = skyinfo.DEFAULT_MINISKY_FILES
            in_filenames = []
            for i in ccdnumlist:
#               in_filenames.append(miniskyfiles.format(i))
                in_filenames.append(miniskyfiles % i)
        
        out_filename = config.get(cls.step_name, 'outfilename')
        logger.info('Sky combine output to %s' % out_filename)
    
        ret_code = cls.__call__(in_filenames, out_filename, mask_value, invalid)
        return ret_code

sky_combine = SkyCombine()

# internal functions & classes

if __name__ == '__main__':
    sky_combine.main()
