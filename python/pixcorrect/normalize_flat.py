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
    def __call__(cls, inlist, ccdnorm):
        """Apply a flat field correction to an image

        :Parameters:
            - `inlist`: list of input and output flat DESImage(s) to normalize
            - `flat_im`:  the flat correction image to apply

        Applies the correction "in place"
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
                        print("Image {:s} did not have a SCALMEAN keyword. Aborting!".format(tmp_dict['fname']))
                        raise
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
            (type, value, trback)=sys.exc_info()
            print("{:s} {:s} {:s} \n".format(inlist,type,value))
            raise
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
                print("CCDNUM: {:d} not among images read in list. Aborting.".format(ccdnorm))
                raise
            logger.info('Normaliztion set to CCD %d value: %.2f ' % (ccdnorm,normval))
        else:
            a_scalmean=np.array(scalmean_list)
            normval=np.median(a_scalmean)
            logger.info('Normaliztion set to median value: %.2f ' % normval )
#
#       Go ahead and normalize the set
#
        logger.info('Normalizing list')
        for tmp_record in norm_list:
            image=DESImage.load(tmp_record['fname'])
            nfactor=tmp_record['scalmean']/normval
            nfactor2=nfactor*nfactor
            logger.info(' CCD: %2d, normalization factor: %.5f ' % (tmp_record['ccdnum'],nfactor) )
            image.data*=nfactor            
            image.weight*=nfactor2            

            ampsecan = [x-1 for x in scan_fits_section(image, 'AMPSECA')]
            ampsecbn = [x-1 for x in scan_fits_section(image, 'AMPSECB')]
            if ((ampsecan[2]>ampsecan[3])or(ampsecbn[2]>ampsecbn[3])):
                print("AMPSEC keyword detected with reversed order")
            print "AMPSECA: ",ampsecan
            print "AMPSECB: ",ampsecbn
#           order x-ranges
            if (ampsecan[0]>ampsecan[1]):
                xtemp=ampsecan[0]
                ampsecan[0]=ampsecan[1]
                ampsecan[1]=xtemp
            if (ampsecbn[0]>ampsecbn[1]):
                xtemp=ampsecbn[0]
                ampsecbn[0]=ampsecbn[1]
                ampsecbn[1]=xtemp
#           order y-ranges
            if (ampsecan[2]>ampsecan[3]):
                ytemp=ampsecan[2]
                ampsecan[2]=ampsecan[3]
                ampsecan[3]=ytemp
            if (ampsecbn[2]>ampsecbn[3]):
                ytemp=ampsecbn[2]
                ampsecbn[2]=ampsecbn[3]
                ampsecbn[3]=ytemp
            print "AMPSECA(ordered): ",ampsecan
            print "AMPSECB(ordered): ",ampsecbn
#

            image.data[ampsecan[2]:ampsecan[3]+1,ampsecan[0]:ampsecan[1]+1]*=image['GAINA']
            image.data[ampsecbn[2]:ampsecbn[3]+1,ampsecbn[0]:ampsecbn[1]+1]*=image['GAINB']

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
#        logger.info('Reading flat correction from %s'% flat_fname)
#        flat_im = DESImage.load(flat_fname)
    
        ret_code = cls.__call__(flat_inlist, ccdnorm)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the flat field correction
        """
        parser.add_argument('--inlist', nargs=1, default=None,
                            help='List of input/output flat field image')
        parser.add_argument('--ccdnorm', type=int, default=-1,
                            help='Specific CCD to use for normalization (default=-1 --> use median over set)')

normalize_flat = NormalizeFlat()

# internal functions & classes

if __name__ == '__main__':
    normalize_flat.main()
