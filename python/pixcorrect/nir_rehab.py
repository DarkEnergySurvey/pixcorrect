#!/usr/bin/env python3
"""Remediate rehab NIR data to fit a DECam-like infrastructure.
   Data currently handled are those taken with VIRCAM class instruemnts
   with standard VSA pre-processing for micro-step and/or dithering 
"""

import time
import os
import sys
import re

import fitsio
import numpy as np

from pixcorrect.corr_util import logger
#from pixcorrect.PixCorrectDriver import PixCorrectMultistep
#from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect.PixCorrectDriver import PixCorrectStep
from pixcorrect.lightbulb_utils import medclip
import pixcorrect.nircaminfo as nci

from despyastro.CCD_corners import get_DESDM_corners_extent
from despyastro import wcsutil, astrometry
from despymisc.misctime import convert_utc_str_to_nite

import despyfits
from despyfits.maskbits import *
from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage, update_hdr_compression, insert_eupspipe
#from despyfits import updateWCS

from astropy.wcs import WCS
import astropy.units as u

# Which section of the config file to read for this step
config_section = 'nirrehab'


class NirRehab(PixCorrectStep):
    description = "Take NIR PawPrint image and confidence files and form a DESImage like set of files"
    step_name = config_section

    @classmethod
    def __call__(cls, nir_paw_image_fname, nir_paw_conf_fname, output_template, conf_limit):
        """Constructs DESIMAGE files from NIR pawprint FITS 

        :Parameters:
            - `nir_paw_image_fname`: filename for the NIR science image (multi-HDU FITS)
            - `nir_paw_conf_fname`:  filename for the NIR confidence image (multi-HDU FITS) 

        """

#       on with the show
        logger.info('Opening science and confidence frames')
        ifits=fitsio.FITS(nir_paw_image_fname,'r')
        cfits=fitsio.FITS(nir_paw_conf_fname,'r')

#
#       Check that the number of HDUs match
#

        if (len(ifits) != len(cfits)):
            print("Number of HDUs/extensions in IMAGE and CONFidence files do not match.")
            print("Aborting")
            exit(1)

        p_ih=ifits[0].read_header()
        p_ch=cfits[0].read_header()
#       Remove reserve keywords
        p_ih.clean()

#
#       Extract some keywords from PRIMARY header to propagate into the individual images.
#
        base_dict={}
        base_header=[]
        for hkeep in nci.nir_paw_primary_keep:
            if (hkeep in p_ih):
                base_header.append({'name':hkeep,'value':p_ih[hkeep],'comment':p_ih.get_comment(hkeep)})
                base_dict[hkeep]={'value':p_ih[hkeep],'comment':p_ih.get_comment(hkeep)}
            else:
                print("Keyword {:s} missing in HDU[{:d}]".format(hkeep,0))
#
#       If possible, need too keep track of REQTIME (requested frametime) because sometimes 
#       EXPTIME seems to be mispopulated in the CCD image HDUs with TEXPTIME
#
        if ('TEXPTIME' in p_ih):
            texptime=p_ih['TEXPTIME']
        else:
            texptime=None
        if ('REQTIME' in p_ih):
            reqtime=p_ih['REQTIME']
        else:
            reqtime=None
#
#        print(base_header)
        

#
#       Step through HDUs... and form "CCD" images for each HDU
#
        ExtList=[]
        for hnum in range(1,len(ifits)):
            print("############ Begin work on extnum={:d} ###############".format(hnum))

#           Check that extensions match (after that concentrate on image).
            print(hnum,ifits[hnum].get_extname(),cfits[hnum].get_extname())
            if (ifits[hnum].get_extname() != cfits[hnum].get_extname()):
                print("Working on extension {:d}.  Extension names (image,conf) of ([{:s}],[{:s}]) do not match!".format(
                    hnum,ifits[hnum].get_extname(),cfits[hnum].get_extname()))
                print("Aborting!")
                exit(1)

            f_ih=ifits[hnum].read_header()
            f_ih.clean()
#
#           Fix occurences where the CCD-level keyword EXPTIME has inherited the value of TEXPTIME
#
            exptime=f_ih['EXPTIME']
            if (reqtime is not None):
                if (exptime > reqtime):
                    print("Warning: possible corrupt EXPTIME (total exptime rather than frame time present).")
                    print("Attempting to update EXTIME to REQTIME (requested frame time).")
                    print("    Primary HDU: TEXPTIME: {:}".format(texptime))
                    print("    Primary HDU: REQTIME:  {:}".format(reqtime))
                    print("    Current HDU: EXPTIME:  {:} --> {:}".format(exptime,reqtime))
                    exptime=reqtime
                    f_ih['EXPTIME']=reqtime
#
#           Augment keywords pulled from primary header with keywords from current HDU
#
            c_header=base_header[:]
            c_dict=dict(base_dict)
            for hkeep in nci.nir_paw_hdu_keep:
                if (hkeep in f_ih):
#                   print(hkeep,f_ih[hkeep],f_ih.get_comment(hkeep))
                    c_header.append({'name':hkeep,'value':f_ih[hkeep],'comment':f_ih.get_comment(hkeep)})
                    if (hkeep in c_dict):
                        print("Warning: Replacing keyword {:s} with value from hdu={:d}".format(hkeep,hnum))
                    c_dict[hkeep]={'value':f_ih[hkeep],'comment':f_ih.get_comment(hkeep)}
                else:
                    print("Keyword {:s} missing in HDU[{:d}]".format(hkeep,hnum))

#
#           Get the CCDNUM from special keyword and propagate
#           Get SKYLEVEL, SKYNOISE, ZEROPOINT and form basis value for the weight plane
#
            ccdnum=f_ih['HIERARCH ESO DET CHIP NO']
            c_header.append({'name':'CCDNUM','value':ccdnum,'comment':'Unique Detector Number'})

#            exptime=f_ih['EXPTIME']
##           Fix occurences where the CCD-level keyword EXPTIME has inherited the value of TEXPTIME
#            if (exptime > reqtime):
#                print("Warning: possible corrupt EXPTIME (total exptime rather than frame time present).")
#                print("Attempting to update EXTIME to REQTIME (requested frame time).")
#                print("    Primary HDU: TEXPTIME: {:.2f}".format(texptime))
#                print("    Primary HDU: REQTIME:  {:.2f}".format(reqtime))
#                print("    Current HDU: EXPTIME:  {:.2f} --> {:.2f}".format(exptime,reqtime))
#                exptime=reqtime
#                f_ih['EXPTIME']=reqtime

            mtime=2.5*np.log10(exptime)
            skylev=f_ih['SKYLEVEL']
            skyrms=f_ih['SKYNOISE']
            seeing=f_ih['SEEING']
            magzpt=f_ih['MAGZPT']

#           zeropoint include a correction from VEGA->AB
#           zeropoint in headers was found to have a factor for EXPTIME removed (have to add back in for DES-like processing)

            if (p_ih['BAND'] in nci.nir_vega_to_ab):
                magzpt=magzpt+nci.nir_vega_to_ab[p_ih['BAND']]+mtime
            else:
                print("Warning! Unknown BAND ({:s}) for conversion of zeropoint from VEGA to AB system".format(p_ih['BAND']))

            c_header.append({'name':'SKYBRITE', 'value':skylev, 'comment':'Sky level estimate from IMCORE'})
            c_header.append({'name':'SKYSIGMA', 'value':skyrms, 'comment':'Sky noise estimate from IMCORE'})
            c_header.append({'name':'SKYVARA',  'value':skyrms*skyrms, 'comment':'Sky noise estimate from IMCORE'})
            c_header.append({'name':'SKYVARB',  'value':skyrms*skyrms, 'comment':'Sky noise estimate from IMCORE'})
            c_header.append({'name':'FWHM',     'value':seeing, 'comment':'Average FWHM (pixels)'})
            c_header.append({'name':'MAG_ZERO', 'value':magzpt, 'comment':'Converted MAGZPT(Vega) to AB system'})
            nite_val=convert_utc_str_to_nite(f_ih['DATE-OBS'])
            c_header.append({'name':'NITE',     'value':nite_val, 'comment':'Observation Nite'})
            c_header.append({'name':'SATURATE', 'value':nci.nircam_satval[ccdnum], 'comment': 'Saturation Level (ADU)'})
            c_header.append({'name':'PIXSCAL1', 'value':0.341, 'comment': 'Fiducial pixel scale (arcsec/pix)'})
            c_header.append({'name':'PIXSCAL2', 'value':0.341, 'comment': 'Fiducial pixel scale (arcsec/pix)'})

#            bval=f_ih['BSCALE']
#            print("BSCALE was: ",bval)
            print("SKYLEVEL was: ",skylev)
            print("SKYRMS was: ",skyrms)
#
#           Searching for a proper WGT prescription
#
#           This was what I took to be equivalent to DES (but perhaps it does not properly factor in N-image stack
#            wgtval=skylev+(skyrms*skyrms)
            print("SKYLEV + (SKYRMS*SKYRMS): ",skylev+(skyrms*skyrms))
#
#           This was assuming SKYLEVEL does not properly inform stats
#            wgtval=(skyrms*skyrms)
            print("(SKYRMS*SKYRMS): ",skyrms*skyrms)

#
#           Read the image data from the science and confidence files.
#
            sci_data=ifits[hnum].read()
            print("Median of data {:.3f} ".format(np.median(sci_data)))
            conf_data=cfits[hnum].read()

#
#           Better seemed to be a re-measurement of STD
#
            print("Attempting an improved SKYRMS with 3-sigma clip to remove objects")
            avgval, medval, stdval = medclip(sci_data,verbose=3)
#            print(avgval,medval,stdval)
            print("stdval^2: ",stdval*stdval)
            wgtval=(stdval*stdval)
#            print(wgtval)
#
#           Use the new (i.e. chip-based header) to feed a WCS 
#           Use image size to feed calculations for center and corners (similar to despyastro.CCD_corners
#
            print("Calculating center/corners assuuming native ZPN projection")
            w=WCS(fitsio.FITSHDR(c_header))

            fnax2=float(sci_data.shape[0])
            fnax1=float(sci_data.shape[1])
            corn_x=np.array([fnax1/2.0,1.,fnax1,fnax1,1.])
            corn_y=np.array([fnax2/2.0,1.,1.,fnax2,fnax2])
            sky = w.pixel_to_world(corn_x,corn_y)
            corn_ra=sky.ra.degree
            corn_dec=sky.dec.degree

            c_header.append({'name':'RA_CENT', 'value':corn_ra[0], 'comment':'RA center'})
            c_header.append({'name':'DEC_CENT','value':corn_dec[0],'comment':'DEC center'})
            for i in range(1,5):
                c_header.append({'name':'RAC{:d}'.format(i), 'value':corn_ra[i], 'comment':'RA corner {:d}'.format(i)})
                c_header.append({'name':'DECC{:d}'.format(i),'value':corn_dec[i],'comment':'DEC corner {:d}'.format(i)})
            RACMIN, RACMAX, DECCMIN, DECCMAX, CROSSRA0 = get_DESDM_corners_extent(corn_ra, corn_dec)
            c_header.append({'name':'RACMIN',  'value':RACMIN,  'comment':'Minimum extent of image in RA'})
            c_header.append({'name':'RACMAX',  'value':RACMAX,  'comment':'Maximum extent of image in RA'})
            c_header.append({'name':'DECCMIN', 'value':DECCMIN, 'comment':'Minimum extent of image in Declination'})
            c_header.append({'name':'DECCMAX', 'value':DECCMAX, 'comment':'Maximum extent of image in Declination'})
            c_header.append({'name':'CROSSRA0','value':CROSSRA0,'comment':'Does Image Span RA 0h (Y/N)'})
            c_header.append({'name':'DESEPOCH','value':'NIREPOCH','comment':'Default DES epoch definition for including NIR data'})
#
#
#
            print("Stripping ZPN projection from WCS and creating a shift to get a rough TAN")
            recs_to_delete=[] 
            for i, hrec in enumerate(c_header):
                if (hrec['name'] == 'CTYPE1'):
                    c_header[i]['value']='RA---TAN'
                if (hrec['name'] == 'CTYPE2'):
                    c_header[i]['value']='DEC--TAN'

                if (hrec['name'] == 'CRVAL1'):
                    c_header[i]['value']=corn_ra[0]
                if (hrec['name'] == 'CRVAL2'):
                    c_header[i]['value']=corn_dec[0]
                if (hrec['name'] == 'CRPIX1'):
                    c_header[i]['value']=fnax1/2.0
                if (hrec['name'] == 'CRPIX2'):
                    c_header[i]['value']=fnax2/2.0

                if (hrec['name'] in ['PV2_1','PV2_2','PV2_3','PV2_4','PV2_5']):
                    recs_to_delete.append(i)
            if (len(recs_to_delete) > 0):
                for i in sorted(recs_to_delete,reverse=True):
                    x=c_header.pop(i)
                    print("Removing: {:}".format(x))

            whack=WCS(fitsio.FITSHDR(c_header))
            skyhack = whack.pixel_to_world(corn_x,corn_y)
            whack_corn_ra=skyhack.ra.degree
            whack_corn_dec=skyhack.dec.degree
            for i in range(5):
                cosdec=np.cos(corn_dec[i]*np.pi/180.)
                dra=3600.*(corn_ra[i]-whack_corn_ra[i])*cosdec
                ddec=3600.*(corn_dec[i]-whack_corn_dec[i])
                print(" WCS shift {:d} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} ".format(ccdnum,corn_ra[i],corn_dec[i],whack_corn_ra[i],whack_corn_dec[i],dra,ddec))

#            for i, hrec in enumerate(c_header):
#                print(i,hrec)

#
#           Form the SCI, MSK, and WGT HDUs
#
            im=DESImage(init_data=True,init_mask=True,init_weight=True,shape=sci_data.shape)

            im.data=np.float32(sci_data)
            msk_wsm=np.where(conf_data<conf_limit)
            im.mask[msk_wsm] |= BADPIX_BPM
            im.weight=np.float32(conf_data/100./wgtval)
#
#           Check for extra conditions where further masking is needed
#           Here is where CCD=6 check was started (now removed and placed 
#               in nir_starmask to take advantage of bright object masking
#


#
#           Deal with individual header-isms and write out SCI, MSK, WGT
#           Note this is using fitsio (duplicating some of the DESIMAGE.save 
#           but customization was needed to deal with foibles of the current
#
            fname=re.sub('%02d','{:02d}'.format(ccdnum),output_template,1)
            ofits = fitsio.FITS(fname, 'rw', clobber=True)

            im.header=fitsio.FITSHDR(c_header)            
            im.header['DES_EXT']='IMAGE'
            im.header = update_hdr_compression(im.header, 'SCI')
            ofits.write(im.data,extname='SCI',header=im.header)


            im.mask_hdr=fitsio.FITSHDR(c_header)            
            im.mask_hdr['DES_EXT']='MASK'
            im.mask_hdr = update_hdr_compression(im.mask_hdr, 'MSK')
            im.mask_hdr['DES_EXT']='MASK'
            ofits.write(im.mask,extname='MSK',header=im.mask_hdr)

#            im.weight_hdr=fitsio.FITSHDR(c_header)            
#            print(im.weight_hdr)
            im.weight_hdr = update_hdr_compression(im.weight_hdr, 'WGT')
#            print(im.weight_hdr)
            im.weight_hdr['DES_EXT']='WEIGHT'
            ofits.write(im.weight,extname='WGT',header=im.weight_hdr)

            ofits.close()
            print("Wrote {:s}".format(fname))
            print(" ")
        

        ifits.close()
        cfits.close()

        ret_code = 0
        return ret_code


    @classmethod
#    def step_run(cls, nir_paw_image, nir_paw_conf, config):
    def step_run(cls, config):
        """Customized execution for taking difference between an image and a comparison

        :Parameters:
            - `nir_paw_image`: multi-HDU FITS file with paw_print science image data
            - `nir_paw_conf`:  multi-HDU FITS file with paw_print confidence image data

        """

        nir_paw_image_fname = config.get(cls.step_name, 'nir_paw_image')
        nir_paw_conf_fname  = config.get(cls.step_name, 'nir_paw_conf')
        output_template     = config.get(cls.step_name, 'output_template')
        conf_limit          = config.getfloat(cls.step_name, 'conf_limit')

#        logger.info('reading Comparison image from %s', comp_fname)
#        comp_im = DESImage.load(comp_fname)

        ret_code = cls.__call__(nir_paw_image_fname, nir_paw_conf_fname, output_template, conf_limit)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to comparison subtraction
        """
        parser.add_argument('--nir_paw_image', nargs=1, default=None,
                            help='NIR pawprint FITS file')
        parser.add_argument('--nir_paw_conf', nargs=1, default=None,
                            help='NIR pawprinti confidence FITS file')
        parser.add_argument('--output_template', nargs=1, default=None,
                            help='Template output filename (%%02d) indicates wildcard for ccdnum/chipnum replacement')
        parser.add_argument('--conf_limit', nargs=1, default=50.,
                            help='Limit in confidence image used to mask low confidece data (default=50.0)')

nir_rehab = NirRehab()

# internal functions & classes

if __name__ == '__main__':
    nir_rehab.main()
