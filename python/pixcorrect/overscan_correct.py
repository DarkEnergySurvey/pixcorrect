#!/usr/bin/env python
"""Apply a overscan correction to a raw DES image 
"""

import ctypes
from os import path
import numpy as np
from pixcorrect import proddir
from pixcorrect.corr_util import logger, load_shlib
from despyfits.DESImage import DESImage, DESImageCStruct, scan_fits_section, data_dtype
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# Which section of the config file to read for this step
config_section = 'overscan'

# Lowest level access to the C library function
#bpm_lib = load_shlib('libbpm')
#bpm_c = bpm_lib.bpm
#bpm_c.restype = ctypes.c_int
#overscan_c.argtypes = [DESImageCStruct, DESImageCStruct]

class ApplyOverScan(PixCorrectImStep):
    description = "Determine and apply an overscan correction"
    step_name = config_section

    @classmethod
    def __call__(cls, image, overscan_sample, overscan_function, overscan_order, overscan_trim):
        """Apply an overscan correction to a raw DES image

        :Parameters:
            - `image`: the DESImage to determine and apply an ovescan correction

        Applies the correction "in place"
        """
 
        logger.info('Applying Overscan')
#        ret_code = bpm_c(image.cstruct, bpm_im.cstruct)
        biassecan = [x-1 for x in scan_fits_section(image, 'BIASSECA')]
        biassecbn = [x-1 for x in scan_fits_section(image, 'BIASSECB')]
        print "BIASSECA: ",biassecan
        print "BIASSECB: ",biassecbn
        biassecan[0]+=overscan_trim         
        biassecan[1]-=overscan_trim         
        biassecbn[0]+=overscan_trim         
        biassecbn[1]-=overscan_trim         
        print "BIASSECA (trimmed): ",biassecan
        print "BIASSECB (trimmed): ",biassecbn

        datasecan = [x-1 for x in scan_fits_section(image, 'DATASECA')]
        datasecbn = [x-1 for x in scan_fits_section(image, 'DATASECB')]
        datasecn = [x-1 for x in scan_fits_section(image, 'DATASEC')]
        print "DATASECA: ",datasecan
        print "DATASECB: ",datasecbn
        if (datasecan[2]>datasecan[3]):
            ytemp=datasecan[2]
            datasecan[2]=datasecan[3]
            datasecan[3]=ytemp
        if (datasecbn[2]>datasecbn[3]):
            ytemp=datasecbn[2]
            datasecbn[2]=datasecbn[3]
            datasecbn[3]=ytemp
        print "DATASECA (ordered): ",datasecan
        print "DATASECB (ordered): ",datasecbn
        print "DATASEC  (ordered): ",datasecn

        ampsecan = [x-1 for x in scan_fits_section(image, 'AMPSECA')]
        ampsecbn = [x-1 for x in scan_fits_section(image, 'AMPSECB')]
        print "AMPSECA: ",ampsecan
        print "AMPSECB: ",ampsecbn
#       order x-ranges
        if (ampsecan[0]>ampsecan[1]):
            xtemp=ampsecan[0]
            ampsecan[0]=ampsecan[1]
            ampsecan[1]=xtemp
        if (ampsecbn[0]>ampsecbn[1]):
            xtemp=ampsecbn[0]
            ampsecbn[0]=ampsecbn[1]
            ampsecbn[1]=xtemp
#       order y-ranges
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
#       Obtain/sample the overscan strip(s) to obtain a set of values.
#
        if (overscan_sample < 0):
            overscan_a=np.median(image.data[biassecan[2]:biassecan[3]+1,biassecan[0]:biassecan[1]+1],axis=1)
            overscan_b=np.median(image.data[biassecbn[2]:biassecbn[3]+1,biassecbn[0]:biassecbn[1]+1],axis=1)
        elif (overscan_sample == 0):
            overscan_a=image.data[biassecan[2]:biassecan[3]+1,biassecan[0]:biassecan[1]+1].mean(axis=1)
            overscan_b=image.data[biassecbn[2]:biassecbn[3]+1,biassecbn[0]:biassecbn[1]+1].mean(axis=1)
        elif (overscan_sample > 0):
            overscan_a=(image.data[biassecan[2]:biassecan[3]+1,biassecan[0]:biassecan[1]+1].sum(axis=1)-
                        image.data[biassecan[2]:biassecan[3]+1,biassecan[0]:biassecan[1]+1].min(axis=1)-
                        image.data[biassecan[2]:biassecan[3]+1,biassecan[0]:biassecan[1]+1].max(axis=1))/float(biassecan[1]-biassecan[0]-1)
            overscan_b=(image.data[biassecbn[2]:biassecbn[3]+1,biassecbn[0]:biassecbn[1]+1].sum(axis=1)-
                        image.data[biassecbn[2]:biassecbn[3]+1,biassecbn[0]:biassecbn[1]+1].min(axis=1)-
                        image.data[biassecbn[2]:biassecbn[3]+1,biassecbn[0]:biassecbn[1]+1].max(axis=1))/float(biassecan[1]-biassecan[0]-1)

        print overscan_a.shape
        print overscan_a
        print overscan_a.dtype
#
#       If requested, reform the overscan estimate using a fit.
#
        if (overscan_function < 0):
#
#           Spline interpolation of resampled set chosen
#
            raise NotImplementedError
        elif (overscan_function > 0):
#
#           Polynomial fit (using Legendre polynomials) of a resampled set chosen
#
            raise NotImplementedError
        elif (overscan_function == 0):
            overscan_val_a=overscan_a 
            overscan_val_b=overscan_b
#
#       Subtract Overscan from the Image data
#
        print datasecn[3]-datasecn[2]+1,datasecn[1]-datasecn[0]+1
        newimage=np.empty([datasecn[3]-datasecn[2]+1,datasecn[1]-datasecn[0]+1],dtype=data_dtype)
        newimage[:,:]=image.data[datasecn[2]:datasecn[3]+1,datasecn[0]:datasecn[1]+1]
        newimage[ampsecan[2]:ampsecan[3]+1,ampsecan[0]:ampsecan[1]+1]-=overscan_val_a[:,np.newaxis]   
        newimage[ampsecbn[2]:ampsecbn[3]+1,ampsecbn[0]:ampsecbn[1]+1]-=overscan_val_b[:,np.newaxis]  
        image.data=newimage
        print image.data.dtype

#
#       Trim data to form output
#
#        print image.data.shape
#        print output.dtype
#        import fitsio
#        fitsio.write('outimage',output,clobber=True)
#        print image.data.shape
#        image.header['NAXIS1']=datasecn[1]-datasecn[0]+1
#        image.header['NAXIS2']=datasecn[3]-datasecn[2]+1

        logger.info(' %d %d ' % image.data.shape)
        logger.debug('Finished applying Overscan')
        ret_code=0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
        """Customized execution for application of the BPM

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
#        import pdb; pdb.set_trace()
#        overscan_apply = config.getboolean(cls.step_name, 'overscan')
        overscan_sample = config.getint(cls.step_name, 'overscansample')
        overscan_function = config.getint(cls.step_name, 'overscanfunction')
        overscan_order = config.getint(cls.step_name, 'overscanorder')
        overscan_trim  = config.getint(cls.step_name, 'overscantrim')
        logger.info('Ovescan will be applied to %s' % image)
    
        ret_code = cls.__call__(image, overscan_sample, overscan_function, overscan_order, overscan_trim)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the BPM
        """
        parser.add_argument('--overscansample', type=int, default=0,
                            help='Overscan sample option (<=-1,0,>=1 MEDIAN, MEAN, MEAN w/minmax reject)')
        parser.add_argument('--overscanfunction', type=int, default=0,
                            help='Overscan function option (<=-1,0,>=1 SPLINE, LINE-by-LINE, Legendre Polynomial)')
        parser.add_argument('--overscanorder', type=int, default=0,
                            help='Order for Legendre polynomial fitting)')
        parser.add_argument('--overscantrim', type=int, default=0,
                            help='Number of ovescan columns to ignore)')


apply_overscan = ApplyOverScan()

# internal functions & classes

if __name__ == '__main__':
    apply_overscan.main()
