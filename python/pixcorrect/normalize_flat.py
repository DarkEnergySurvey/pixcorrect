#!/usr/bin/env python
"""Apply a flat correction to a raw DES image 
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
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'normflat'

# Lowest level access to the C library function
#flatcorrect_lib = load_shlib('libflatcorrect')
#flat_c = flatcorrect_lib.flat_c
#flat_c.restype = ctypes.c_int
#flat_c.argtypes = [DESImageCStruct, DESImageCStruct]

class NormalizeFlat(PixCorrectImStep):
    description = "Normalize a set of flat fields"
    step_name = config_section

    @classmethod
    def __call__(cls, inlist, ccdnorm, ampborder):
        """Apply a flat field correction to an image

        :Parameters:
            - `inlist`: list of input and output flat DESImage(s) to normalize
            - `flat_im`:  the flat correction image to apply

        Applies the correction to each input and writes a separate output file.
        """
 
        logger.info('Initial Read of Flat Field Headers')
#
        norm_list=[]
        scalmean_list=[]
        normval=None
#
        try:
            f1=open(inlist,'r')
            for line in f1:
                line=line.strip()
                columns=line.split()
                if (os.path.isfile(columns[0])):    
                    tmp_dict={}
                    tmp_dict['fname']=columns[0]
                    tmp_dict['oname']=columns[1]
                    if (tmp_dict['fname'][-2:] == "fz"):
                        sci_hdu=1 # for .fz
                    else:
                        sci_hdu=0 # for .fits (or .gz)
                    temp_fits=fitsio.FITS(tmp_dict['fname'],'r')
                    temp_head=temp_fits[sci_hdu].read_header()
#
#                   Get the CCD number
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
#                   Get the SCALMEAN value
#
                    try:
                        tmp_dict['scalmean']=float(temp_head['SCALMEAN'])
                    except:
                        raise ValueError("Image %s did not have a SCALMEAN keyword. Aborting!" % tmp_dict['fname'])
#
#                   Finished first header census
#                   Save file info and scalmean's to a list
#
                    norm_list.append(tmp_dict)
                    scalmean_list.append(tmp_dict['scalmean'])
                    temp_fits.close()
            f1.close()
        except:
#
#           Input file was not present.           
#
#            (type, value, trback)=sys.exc_info()
#            print("{:s} {:s} {:s} \n".format(inlist,type,value))
            raise IOError("File not found.  Missing input list %s " % inlist )
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
#       Go ahead and normalize the set
#
        logger.info('Normalizing list')
        for tmp_record in norm_list:
            logger.info('Working on image: %s ' % (tmp_record['fname']) )
            image=DESImage.load(tmp_record['fname'])
            nfactor=tmp_record['scalmean']/normval
            nfactor2=nfactor*nfactor
            logger.info(' CCD: %2d, relative normalization factor: %.5f ' % (tmp_record['ccdnum'],nfactor) )
            image.data*=nfactor            
            image.weight*=nfactor2            
#
#           Create keywords that reflect the median value of the flat on each amp.
#
            for amp in decaminfo.amps:
                datasecn=scan_fits_section(image,'DATASEC'+amp)
                datasecn[0]=datasecn[0]+ampborder
                datasecn[1]=datasecn[1]-ampborder
                datasecn[2]=datasecn[2]+ampborder
                datasecn[3]=datasecn[3]-ampborder
                image['FLATMED'+amp]=np.median(image.data[datasecn[2]:datasecn[3]+1,datasecn[0]:datasecn[1]+1])
            
            DESImage.save(image,tmp_record['oname'])

        logger.debug('Finished applying Flat')
        ret_code=0

        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the Bias

        :Parameters:
            - `image`: the DESImage on which to operate
            - `flat`: the bias image to apply

        """

        flat_inlist = config.get(cls.step_name, 'inlist')
        ccdnorm = config.getint(cls.step_name, 'ccdnorm')
        ampborder = config.getint(cls.step_name, 'ampborder')

        ret_code = cls.__call__(flat_inlist, ccdnorm, ampborder)

        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--inlist', nargs=1, default=None,
                            help='List of input/output flat field image')
        parser.add_argument('--ccdnorm', type=int, default=-1,
                            help='Specific CCD to use for normalization (default=-1 --> use median over set)')
        parser.add_argument('--ampborder', type=int, default=50,
                            help='Length in pixels around periphery of each amp to ignore when calculating statistics (default=50)')

normalize_flat = NormalizeFlat()

# internal functions & classes

if __name__ == '__main__':
    normalize_flat.main()
