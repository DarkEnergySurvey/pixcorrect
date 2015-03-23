#!/usr/bin/env python
"""Find a normalization for a set of flat field images
"""

import ctypes
import sys
import os
#from os import path
import fitsio
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectStep, filelist_to_list
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'findflatnorm'

# Lowest level access to the C library function
#flatcorrect_lib = load_shlib('libflatcorrect')
#flat_c = flatcorrect_lib.flat_c
#flat_c.restype = ctypes.c_int
#flat_c.argtypes = [DESImageCStruct, DESImageCStruct]

class FindFlatNormalization(PixCorrectStep):
    description = "Find normalization factor for a set of flat fields"
    step_name = config_section

    @classmethod
    def __call__(cls, in_filenames, ccdnorm, outnorm):
        """Apply a flat field correction to an image

        :Parameters:
            - `in_filenames`: list of input DESImage(s) to use to determine the normalization factor 
            - `ccdnorm`: -1-->normalize to median of all files, or to image with CCDNUM=ccdnorm 
            - `outnorm`: output file name to write the normalization factor

        Applies the correction to each input and writes a separate output file.
        """
 
        logger.info('Initial Read of Flat Field Headers')
#
        norm_list=[]
        scalmean_list=[]
        normval=None
#
        for filename in in_filenames:
            if (os.path.isfile(filename)):    
                tmp_dict={}
                tmp_dict['fname']=filename
                if (tmp_dict['fname'][-2:] == "fz"):
                    sci_hdu=1 # for .fz
                else:
                    sci_hdu=0 # for .fits (or .gz)
                temp_fits=fitsio.FITS(tmp_dict['fname'],'r')
                temp_head=temp_fits[sci_hdu].read_header()
#
#               Get the CCD number
#
                try:
                    tmp_dict['ccdnum']=int(temp_head['CCDNUM'])
                except:
                    if (ccdnorm < 1):
                        tmp_dict['ccdnum']=-1
                        pass
                    else:
                        print("Warning: image {:s} did not have a CCDNUM keyword!".format(tmp_dict['fname']))
                        pass
#
#               Get the SCALMEAN value
#
                try:
                    tmp_dict['scalmean']=float(temp_head['SCALMEAN'])
                except:
                    raise ValueError("Image %s did not have a SCALMEAN keyword. Aborting!" % tmp_dict['fname'])
#
#               Finished first header census
#               Save file info and scalmean's to a list
#
                norm_list.append(tmp_dict)
                scalmean_list.append(tmp_dict['scalmean'])
                temp_fits.close()
#
#       All information is now present. Determine the value that will be used in normalization.
#
        if (ccdnorm > 1):
            for tmp_rec in norm_list:
                if (normval is None):
                    if (tmp_rec['ccdnum']==ccdnorm):
                        normval=tmp_rec['ccdnum']
                else:
                    if (tmp_rec['ccdnum']==ccdnorm):
                        print("Warning: More than one image with CCDNUM={:d} identified")
            if (normval is None):
                raise ValueError("No image with CCDNUM=%d found among input list. Aborting!" % ccdnorm)
            logger.info('Normaliztion: %.2f set based on value from CCD %d ' % (normval,ccdnorm))
        else:
            a_scalmean=np.array(scalmean_list)
            normval=np.median(a_scalmean)
            logger.info('Normaliztion: %.2f set based on median value of the ensemble ' % normval )
#
#       Write out the normalization factor
#
        fout=open(outnorm,'w')
        fout.write("{:.2f}\n".format(normval))
        fout.close()
        ret_code=0

        return ret_code


    @classmethod
    def step_run(cls, config):
        """Customized execution to find normalization factor

        :Parameters:
            - `image`: the DESImage on which to operate

        """

        flat_inlist = config.get(cls.step_name, 'inlist')
        in_filenames=filelist_to_list(flat_inlist)
        ccdnorm = config.getint(cls.step_name, 'ccdnorm')
        outnorm = config.get(cls.step_name, 'outnorm')

        ret_code = cls.__call__(in_filenames, ccdnorm, outnorm)

        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--inlist', nargs=1, default=None,
                            help='List of input/output flat field image')
        parser.add_argument('--ccdnorm', type=int, default=-1,
                            help='Specific CCD to use for normalization (default=-1 --> use median over set)')
        parser.add_argument('--outnorm', nargs=1, default=None,
                            help='Output file to write normalization')

find_flat_normalization = FindFlatNormalization()

# internal functions & classes

if __name__ == '__main__':
    find_flat_normalization.main()
