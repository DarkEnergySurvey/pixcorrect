#!/usr/bin/env python3

# $Id: cti.py 47952 2019-01-03 21:04:53Z rgruendl $
# $Rev:: 47952                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2019-01-03 15:04:53 #$:  # Date of last commit.

"""Perform Masking based on a catalog
"""

import os
import numpy as np

import despydb.desdbi
#from despyastro import wcsutil, astrometry
from despyastro import astrometry
from esutil     import wcsutil

from pixcorrect import starmask_util as smu
from pixcorrect.lightbulb_utils import medclip
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from despyfits.DESImage import DESImage, DESImageCStruct, section2slice, data_dtype
import pixcorrect.nircaminfo as nci
from pixcorrect import region_utils as ru
from pixcorrect.clip_mask_utils import polygon_to_pix

import despyfits
from despyfits.maskbits import *
from despyfits.maskbits import parse_badpix_mask


# Which section of the config file to read for this step
config_section = 'nir_starmask'

class NIRStarMask(PixCorrectImStep):
    description = "Mask stars in NIR image based on a recipe"
    step_name = config_section

    @classmethod
    def __call__(cls, image, dbSection, UseBand, SE_Special, RegionFile, RegionFlagVal):
        """
        Read a FITS catalog of 2MASS PSC sources and mask areas around bright sources based
        on their magnitude.
        """

        verbose=2

        # Populate radec_box structure with image extent
        # Expand box with the size of a mask for a zero mag star (~4 arcmin)
        
        radec_box={}
        if (image.header['CROSSRA0']=="Y"):
            radec_box['crossra0']=True
        else:
            radec_box['crossra0']=False
        radec_box['ra1']=image.header['RACMIN']
        radec_box['ra2']=image.header['RACMAX']
        radec_box['dec1']=image.header['DECCMIN']
        radec_box['dec2']=image.header['DECCMAX']

        radec_box=smu.expand_range(radec_box,method='fixed',extend=4.,verbose=verbose)
#
#       Welcome to the pixel-scale rabbit hole.
#       -1- First try to simply get the pixel scale based upon the DES (single-epoch) PIXSCAL1/2 keywords
#           -2- Failing that look for old-school CDELT definitions
#               -3- Failing futher go to the definition based on a CD-matrix
#                   -4- Then the oddball CD-matrix that thinks because its off-diaganol terms are zero it can just omit them
#                   
        try:
            pix_scale=0.5*(image.header['PIXSCAL1']+image.header['PIXSCAL2'])
        except KeyError:
            print("Missing PIXSCAL1/2 keywords?  Attempting with CDELT1/2.")
            try:
                pix_scale=3600.*0.5*(np.abs(image.header['CDELT1'])+np.abs(image.header['CDELT2']))
            except KeyError:
                print("No CDELT1/2 keywords.  Attempting with CD-Matrix")
                try:
                    pix_scale=3600.*0.5*(np.sqrt((image.header['CD1_1']*image.header['CD1_1'])+(image.header['CD2_1']*image.header['CD2_1']))+
                                         np.sqrt((image.header['CD1_2']*image.header['CD1_2'])+(image.header['CD2_2']*image.header['CD2_2'])))
                except KeyError:
                    print("Perhaps not a full CD-Matrix? checking for diagonal elements")
                    try:
                        pix_scale=3600.*0.5*(np.abs(image.header['CD1_1'])+np.abs(image.header['CD2_2']))
                    except Exception as e:
                        print("Exception {:}".format(e))
                        print("Could not determine a pixel scale!  Aborting!")
                        exit(1)
        print("Pixel Scale being used is: {:} arcsec/pixel ".format(pix_scale))
#
#       Get the WCS
#
        wcs=wcsutil.WCS(image.header)

#
#       Prepare a database connection to obtain catalog info
#       Make an RA-Dec box query of 2MASS and get positions with J-band magnitudes.
#
        try:
            desdmfile = os.environ["des_services"]
        except KeyError:
            desdmfile = None
        dbh = despydb.desdbi.DesDbi(desdmfile, dbSection, retry=True)

#
#       Where to get the catalog of stars to mask.
#
        useGAIA=False
        use2MASS=False
        if (UseBand):
            if ('BAND' in image.header):
                print("Checked Image for BAND and found: {:s}".format(image.header['BAND']))
                if (image.header['BAND'] in ['VY','J','H','Ks']):
                    use2MASS=True
                else:
                    useGAIA=True
            else:
                use2MASS=True
        else:
            use2MASS=True

        if (use2MASS):
            print("Querying 2MASS_PSC")
            StarCat, StarHead=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='TWOMASS_PSC',cols=['RA','DEC','J'],Timing=False,verbose=0)
            CatMag='J'
        else:
            print("Querying GAIA_DR2")
            StarCat, StarHead=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='GAIA_DR2',cols=['RA','DEC','PHOT_G_MEAN_MAG'],Timing=False,verbose=0)
            CatMag='PHOT_G_MEAN_MAG'

#
#       Establish bounds of regions to be masked.  And magnitude limit where masking will occur...
#
        if (use2MASS):
#
#           For 2MASS (VIRCAM observations):
#               Current cutoff is J=14.8 (roughly a J-band star that is saturated for 30-s exposures)
#               could also use an exptime+band based calculation (an estimate is provided but not used).
#
            mag2maskrad = [ 0.00022263, -0.00726182,  0.06189586]
            rmask=np.polyval(mag2maskrad,StarCat[CatMag])
#
#           For VISTA VY- and J-band images the ghosting is known to be more severe
#           Therefore, in cases where the relation below gives a larger radius... replace with larger radius for J- and VY-band images.
#
            if ('BAND' in image.header):
                if (image.header['BAND'] in ['VY','J']):
                    polf_coeff =np.array([ 8.97832817e-05, -2.60639835e-03, 1.86773994e-02, 6.42125903e-03])
                    rmask2=np.polyval(polf_coeff,StarCat[CatMag])
                    wsm=np.where(rmask2>rmask)
                    rmask[wsm]=rmask2[wsm]

            if (image.header['BAND'] in nci.nircam_satcalc):
                mag_sat_calc=nci.nircam_satcalc[image.header['BAND']]+2.5*np.log10(image.header['EXPTIME'])
            else:
                mag_sat_calc=14.8
#
            mag_sat=14.8 
            print("Using a magnitude cutoff of {:.2f}.  (Note exptime/band based saturation calculation for {:}-band @ {:.1f} sec exposure gives: {:.2f}".format(
                mag_sat,image.header['BAND'],image.header['EXPTIME'],mag_sat_calc))

        else:
#
#           For GAIA (if chosen)
#               Rough Cutoff that is used is GAIA G = 16. (for lack of a better choice)
#               THIS IS ABSOLUTELY A PLACEHOLDER FOR A REAL relationship and cutoff....
#
            print("WARNING:  The code is attempting to use GAIA but the radius relation is merely that used for 2MASS J and typical VIRCAM observing modes.")
            mag2maskrad = [ 0.00022263, -0.00726182,  0.06189586]
#            print(np.polyval(mag2maskrad,0.0))
            rmask=np.polyval(mag2maskrad,StarCat[CatMag])
            mag_sat=16.

#
#       Some code that can be set to dump a DS9 region file if someone wanted to make a quick comparison.
#   
        dump_regions=False 
        if (dump_regions):
            print("##################################")
            print("#   2MASS PSC based mask regions  ")
            for i in range(rmask.size):
                if (StarCat['CatMag'][i]<mag_sat):
                    print(" fk5;circle({:.6f},{:.6f},{:.2f}\") # color=blue width=2 ".format(StarCat['RA'][i],StarCat['DEC'][i],rmask[i]*3600.))
                else:
                    print(" fk5;circle({:.6f},{:.6f},{:.2f}\") # color=red  width=2 ".format(StarCat['RA'][i],StarCat['DEC'][i],rmask[i]*3600.))

        wsm=np.where(StarCat[CatMag]<mag_sat)
        if (rmask[wsm].size > 0):
            x_cen,y_cen=wcs.sky2image(StarCat['RA'][wsm],StarCat['DEC'][wsm])
            r_pix=rmask[wsm]*3600./pix_scale

            ix1=x_cen-r_pix 
            ix2=x_cen+r_pix 
            iy1=y_cen-r_pix 
            iy2=y_cen+r_pix 

            # Switched from NAXIS keyword to prevent confusion when working with compresssed files
            # nx,ny needed to make quick checks to avoid wasting time on objects where the mask won't overlap
            (ny,nx)=image.mask.shape
#            print("NAXIS (shape): {:d},{:d}".format(nx,ny))
#           Form index arrays for rapid slicing to from circular masks.
            y, x = np.indices(image.mask.shape)
            for i in range(ix1.size):
#
#               Check that it is even possible for the star mask to intersect the image before doing the work
#
                if ((ix2[i]<0)or(ix1[i]>nx)or(iy2[i]<0)or(iy1[i]>ny)):
                    if (verbose>0):
                        print("Skipping (ra,dec  (x,y)  mag,rad):       {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f} ".format(
                            StarCat['RA'][wsm][i],StarCat['DEC'][wsm][i],x_cen[i],y_cen[i],StarCat[CatMag][wsm][i],r_pix[i]))

                else:
                    r = np.sqrt((x - x_cen[i])**2 + (y - y_cen[i])**2)
                    r_wsm=np.where(r<r_pix[i])
                    image.mask[y[r_wsm],x[r_wsm]] |= BADPIX_STAR

                    print("Masking  (ra,dec  (x,y)  mag,rad  Npix): {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f}  {:d} ".format(
                        StarCat['RA'][wsm][i],StarCat['DEC'][wsm][i],x_cen[i],y_cen[i],StarCat[CatMag][wsm][i],r_pix[i],r[r_wsm].size))

#
#       Region_Mask chosen
#
        if (RegionFile):
            reg_dict=ru.loadRegionFile(RegionFile)
            (ny,nx)=image.mask.shape
            print("Identified {:d} regions to be associated with this image".format(reg_dict['nentry']))
            nr_flag_tot=0
            for ireg in range(1,reg_dict['nentry']+1):
#                print(reg_dict[ireg]['line'])
                x, y = wcs.sky2image(reg_dict[ireg]['ra'],reg_dict[ireg]['dec'])
#                print(x,y)
                if (reg_dict[ireg]['type'] == "point"):
#                    print("point")
                    xmsk=np.rint(x)+xpt_exp
                    ymsk=np.rint(y)+ypt_exp
                    xy = np.array([(xmsk[i], ymsk[i]) for i in range(x.size)],dtype=int)
                elif (reg_dict[ireg]['type'] == "polygon"):
#                    print("polygon")
                    xy = polygon_to_pix(x,y)
#                print(xy)

#               # Check for and remove any portion of the region that would extend beyond the image boundaries
                wsm=np.logical_and(np.logical_and(xy[:,0]-1>=0,xy[:,0]<nx),np.logical_and(xy[:,1]-1>=0,xy[:,1]<ny))
                # Flag remainging pixels
                image.mask[xy[:,1][wsm]-1,xy[:,0][wsm]-1] |= parse_badpix_mask(RegionFlagVal)
                nr_flag=xy[:,1][wsm].size
                nr_flag_tot=nr_flag_tot+nr_flag
                # print("Region flagged {:d} pixels".format(nr_flag))
            print("Finished flagging based on region file.  Total of {:d} (not necessarily unique) pixels flagged as {:s}.".format(nr_flag_tot,RegionFlagVal))

#
#       Quick check for CCD=6, readout 14 failure for VISTA images
#       https://www.eso.org/observing/dfo/quality/VIRCAM/pipeline/problems.html
#
        if (SE_Special):
#
#           Perform a series of preliminary checks before even trying to look at the image and mask values
#
            print("Checking whether image is a candidate to check for VISTA ccd6, readout 14 failure")
            SE_Special_DoIt=True
            if ('BAND' in image.header):
                if (not(image.header['BAND'] in ['VY','J','H','Ks'])):
                    SE_Special_DoIt=False
                    print("BAND is not a VISTA NIR band.  Skipping...")
                else:
                    print("Verified VISTA NIR band: {:s}".format(image.header['BAND']))
            else:
                SE_Special_DoIt=False
                print("BAND keyword not present.  Skipping...")
#
            if ('CCDNUM' in image.header):
                if (image.header['CCDNUM']!=6):
                    SE_Special_DoIt=False
                    print("CCDNUM is not #6.  Skipping...")
                else:
                    print("Verified CCD={:d}".format(image.header['CCDNUM']))
            else:
                SE_Special_DoIt=False
                print("CCDNUM keyword not present.  Skipping...")
#
            if ('NITE' in image.header):
                nite_val=int(image.header['NITE'])
                if ((nite_val < 20091015)or(nite_val > 20091120)):
                    SE_Special_DoIt=False
                    print("NITE is not between 20091014 and 20091121 (was {:d}).  Skipping...".format(nite_val))
                else:
                    print("Verified NITE={:d}".format(nite_val))
            else:
                SE_Special_DoIt=False
                print("NITE keyword not present.  Skipping...")
#
            if (SE_Special_DoIt):
                print("Preliminary checks passed... checking outlier rate in suspect region.")

                avgval, medval, stdval = medclip(image.data,verbose=3)

                m3=[medval-3.*stdval,medval+3.*stdval]
                m5=[medval-5.*stdval,medval+5.*stdval]

                (ny,nx)=image.data.shape
                smear=nx-2048
                if (smear < 0):
                    smear=0

                ReadOutList=[11,12,14]
                n3o_sz={}
                n5o_sz={}
                for ir in ReadOutList:
                    ix1=int(((ir-1)*128))
                    ix2=int(((ir)*128)+smear)
                    mwsm=np.where(image.mask[:,ix1:ix2]==0)
                    r_sz=image.data[mwsm].size
                    c14wsm=np.where(np.logical_and(image.mask[:,ix1:ix2]==0,np.logical_or(image.data[:,ix1:ix2]<m3[0],image.data[:,ix1:ix2]>m3[1])))
                    n3o_sz[ir]=100.*image.data[:,ix1:ix2][c14wsm].size/r_sz
                    c14wsm=np.where(np.logical_and(image.mask[:,ix1:ix2]==0,np.logical_or(image.data[:,ix1:ix2]<m5[0],image.data[:,ix1:ix2]>m5[1])))
                    n5o_sz[ir]=100.*image.data[:,ix1:ix2][c14wsm].size/r_sz

                maskReadout=False
                if ((n3o_sz[14]>1.5)and(n5o_sz[14]>1.0)):
                    print("  Readout14 %-3sigma-outliers: {:5.2f}  (criterion is >1.5)".format(n3o_sz[14]))
                    print("  Readout14 %-5sigma-outliers: {:5.2f}  (criterion is >1.0)".format(n5o_sz[14]))
                    print("  Checking relative number to those in readout 11 and 12 areas")
                    r11_3o=n3o_sz[14]/n3o_sz[11]
                    r11_5o=n5o_sz[14]/n5o_sz[11]
                    r12_3o=n3o_sz[14]/n3o_sz[12]
                    r12_5o=n5o_sz[14]/n5o_sz[12]
                    for ir in [11,12]:
                        print("  R14/R{:02d} ratio of 3sigma-outliers: {:7.2f}  (criterion is >3.0)".format(ir,n3o_sz[14]/n3o_sz[ir]))
                        print("  R14/R{:02d} ratio of 5sigma-outliers: {:7.2f}  (criterion is >3.0)".format(ir,n5o_sz[14]/n5o_sz[ir]))
                    if ((r11_3o>3.0)or(r11_5o>3.0)or(r12_3o>3.0)or(r12_5o>3.0)):
                        maskReadout=True
                if (maskReadout):
                    print("Masking Readout14 with BADPIX_BADAMP {:d}".format(BADPIX_BADAMP))
                    ir=14
                    ix1=int(((ir-1)*128))
                    ix2=int(((ir)*128)+smear)
                    image.mask[:,ix1:ix2] |= BADPIX_BADAMP

#       End Special Single-Epoch Mask for Readout 14 of CCD=6 of VISTA/NIR frames                    

        ret_code = 0
        return ret_code


    @classmethod
    def step_run(cls, image, config):
#    def step_run(cls, config):
        """Customized execution for NIR star masking

        :Parameters:
            - `image`: the DESImage on which to operate
            - `section`: the DB section in .desservices.ini to obtain connection info
#            - `config`: the configuration from which to get other parameters

        """
#        logger.info('NIR Starmask {:s}'.format(image))

        dbSection = config.get(cls.step_name, 'section')
        useband  = config.getboolean(cls.step_name, 'useband')
        se_special = config.getboolean(cls.step_name, 'se_special')
        region_file =  config.get(cls.step_name, 'region_mask')
        region_flag_val = config.get(cls.step_name, 'reg_flag_val')     

        ret_code = cls.__call__(image,dbSection,useband,se_special,region_file,region_flag_val)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the nir_starmask utility
        """
        parser.add_argument('--section', nargs=1, default='db-dessci',
                            help='section of .desservices file with connection info')
        parser.add_argument('--useband', action='store_true', default=False, 
                            help='use BAND keyword to choose masking catalog and radius relation')
        parser.add_argument('--se_special', action='store_true', default=False, 
                            help='perform special single-epoch masking (e.g. ccd6 readout14 failure')
        parser.add_argument('--region_mask', action='store', type=str, default='',
                            help='Region file defining masks to be applied to image')
        parser.add_argument('--reg_flag_val', action='store', type=str, default="TRAIL",
                            help='DESDM BADPIX flag bit to apply to pixels in region_mask (default=TRAIL)')

#        parser.add_argument('--Schema', nargs=1, default='des_admin',
#                            help='DB schema (do not include \'.\').')


nir_starmask = NIRStarMask()

# internal functions & classes

if __name__ == '__main__':
    nir_starmask.main()
