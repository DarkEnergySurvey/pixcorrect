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
from despyastro import wcsutil, astrometry

from pixcorrect import starmask_util as smu
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from despyfits.DESImage import DESImage, DESImageCStruct, section2slice, data_dtype
import pixcorrect.nircaminfo as nci

import despyfits
from despyfits.maskbits import *
from despyfits.maskbits import parse_badpix_mask


# Which section of the config file to read for this step
config_section = 'nir_starmask'

class NIRStarMask(PixCorrectImStep):
    description = "Mask stars in NIR image based on a recipe"
    step_name = config_section

    @classmethod
    def __call__(cls, image, dbSection, UseBand):
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
#            print(np.polyval(mag2maskrad,0.0))
            rmask=np.polyval(mag2maskrad,StarCat[CatMag])

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
                if ((ix2[i]<0)or(ix1[i]>nx)or(iy2[i]<0)or(iy2[i]>ny)):
                    if (verbose>0):
                        print("Skipping (ra,dec  (x,y)  mag,rad):       {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f} ".format(
                            StarCat['RA'][wsm][i],StarCat['DEC'][wsm][i],x_cen[i],y_cen[i],StarCat[CatMag][wsm][i],r_pix[i]))

                else:
                    r = np.sqrt((x - x_cen[i])**2 + (y - y_cen[i])**2)
                    r_wsm=np.where(r<r_pix[i])
                    image.mask[y[r_wsm],x[r_wsm]] |= BADPIX_STAR

                    print("Masking  (ra,dec  (x,y)  mag,rad  Npix): {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f}  {:d} ".format(
                        StarCat['RA'][wsm][i],StarCat['DEC'][wsm][i],x_cen[i],y_cen[i],StarCat[CatMag][wsm][i],r_pix[i],r[r_wsm].size))

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

        ret_code = cls.__call__(image,dbSection,useband)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """
        parser.add_argument('--section', nargs=1, default='db-dessci',
                            help='section of .desservices file with connection info')
        parser.add_argument('--useband', action='store_true', default=False, 
                            help='use BAND keyword to choose masking catalog and radius relation')

#        parser.add_argument('--Schema', nargs=1, default='des_admin',
#                            help='DB schema (do not include \'.\').')


nir_starmask = NIRStarMask()

# internal functions & classes

if __name__ == '__main__':
    nir_starmask.main()
