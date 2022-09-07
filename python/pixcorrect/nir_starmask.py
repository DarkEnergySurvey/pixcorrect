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
    def __call__(cls, image, dbSection):
        """
        Read a FITS catalog of 2MASS PSC sources and mask areas around bright sources based
        on their magnitude.
        """

        verbose=2

        # Populate radec_box structure with image extent
        # Expand box with the size of a mask for a zero mag star (~4 arcmin)
        
        radec_box={}
        if (image['CROSSRA0']=="Y"):
            radec_box['crossra0']=True
        else:
            radec_box['crossra0']=False
        radec_box['ra1']=image['RACMIN']
        radec_box['ra2']=image['RACMAX']
        radec_box['dec1']=image['DECCMIN']
        radec_box['dec2']=image['DECCMAX']

        radec_box=smu.expand_range(radec_box,method='fixed',extend=4.,verbose=verbose)
        pix_scale=0.5*(image['PIXSCAL1']+image['PIXSCAL2'])
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

        cat_2mass, header_2mass=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='TWOMASS_PSC',cols=['RA','DEC','J'],Timing=False,verbose=0)

#
#       Establish bounds of regions to be masked.
#
        mag2maskrad = [ 0.00022263, -0.00726182,  0.06189586]
        print(np.polyval(mag2maskrad,0.0))
        rmask=np.polyval(mag2maskrad,cat_2mass['J'])

#
#       Current cutoff is J=14.8 (roughly a J-band star that is saturated for 30-s exposures)
#       could also use an exptime+band based calculation
#
        if (image['BAND'] in nci.nircam_satcalc):
            mag_sat_calc=nci.nircam_satcalc[image['BAND']]+2.5*np.log10(image['EXPTIME'])
        else:
            mag_sat_calc=14.8
#
        mag_sat=14.8 
        print("Using a magnitude cutoff of {:.2f}.  (Note exptime/band based saturation calculation for {:}-band @ {:.1f} sec exposure gives: {:.2f}".format(
            mag_sat,image['BAND'],image['EXPTIME'],mag_sat_calc))
   
        dump_regions=False 
        if (dump_regions):
            print("##################################")
            print("#   2MASS PSC based mask regions  ")
            for i in range(rmask.size):
                if (cat_2mass['J'][i]<mag_sat):
                    print(" fk5;circle({:.6f},{:.6f},{:.2f}\") # color=blue width=2 ".format(cat_2mass['RA'][i],cat_2mass['DEC'][i],rmask[i]*3600.))
                else:
                    print(" fk5;circle({:.6f},{:.6f},{:.2f}\") # color=red  width=2 ".format(cat_2mass['RA'][i],cat_2mass['DEC'][i],rmask[i]*3600.))

        wsm=np.where(cat_2mass['J']<mag_sat)
        if (rmask[wsm].size > 0):
            x_cen,y_cen=wcs.sky2image(cat_2mass['RA'][wsm],cat_2mass['DEC'][wsm])
            r_pix=rmask[wsm]*3600./pix_scale

            ix1=x_cen-r_pix 
            ix2=x_cen+r_pix 
            iy1=y_cen-r_pix 
            iy2=y_cen+r_pix 

            nx=image['NAXIS1']
            ny=image['NAXIS2']
            y, x = np.indices(image.mask.shape)
            for i in range(ix1.size):
#
#               Check that it is even possible for the star mask to intersect the image before doing the work
#
                if ((ix2[i]<0)or(ix1[i]>nx)or(iy2[i]<0)or(iy2[i]>ny)):
                    if (verbose>0):
                        print("Skipping (ra,dec  (x,y)  mag,rad):       {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f} ".format(
                            cat_2mass['RA'][wsm][i],cat_2mass['DEC'][wsm][i],x_cen[i],y_cen[i],cat_2mass['J'][wsm][i],r_pix[i]))

                else:
                    r = np.sqrt((x - x_cen[i])**2 + (y - y_cen[i])**2)
                    r_wsm=np.where(r<r_pix[i])
                    image.mask[y[r_wsm],x[r_wsm]] |= BADPIX_STAR

                    print("Masking  (ra,dec  (x,y)  mag,rad  Npix): {:8.5f} {:8.5f}  {:7.1f} {:7.1f}  {:6.2f} {:6.2f}  {:d} ".format(
                        cat_2mass['RA'][wsm][i],cat_2mass['DEC'][wsm][i],x_cen[i],y_cen[i],cat_2mass['J'][wsm][i],r_pix[i],r[r_wsm].size))


#                ctiDict = cti.check_cti(image, CTI[image['CCDNUM']], verbose=1)
##
##               Current criterion:
##                   Looks for horizontal striping in image (with large deficits in counts that are not
##                       associated with an edge-bleed.
##                   Examines auto-correlation for lags in the x-direction at 5, 7, and 15 pixel offsets
##                       and compares to lags obtained from measurments in the diaganol direction.
##                   Looks for evidence of excessive power in the ratio between x-direction and diagnol sets
##                       that indicative that charge is bleeding in the x-direction.
##
#                if ctiDict['isCTI']:
#                    image = cti.mask_cti(image, CTI[image['CCDNUM']], ctiDict, verbose=1)
#                    logger.info(' CTI: Detected CTI for Exp={:d}, CCD={:d}, Amp={:s}'.format(image['EXPNUM'], image['CCDNUM'], CTI[image['CCDNUM']]['amp']))
#                    image.write_key('DES_CTI', 'Masked DATASEC{:s}'.format(CTI[image['CCDNUM']]['amp']))
#
#        logger.debug('Finished checking and applying mask CTI')
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

        ret_code = cls.__call__(image,dbSection)
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific application of the gain correction
        """
        parser.add_argument('--section', nargs=1, default='db-dessci',
                            help='section of .desservices file with connection info')
#        parser.add_argument('--Schema', nargs=1, default='des_admin',
#                            help='DB schema (do not include \'.\').')


nir_starmask = NIRStarMask()

# internal functions & classes

if __name__ == '__main__':
    nir_starmask.main()
