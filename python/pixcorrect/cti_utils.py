#!/usr/bin/env python3

# $Id: cti_utils.py 47952 2019-01-03 21:04:53Z rgruendl $
# $Rev:: 47952                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2019-01-03 15:04:53 #$:  # Date of last commit.

"""CTI masking functions
"""

import numpy as np
from despyfits.DESImage import section2slice
from despyfits.maskbits import *
from pixcorrect import lightbulb_utils as lb
from pixcorrect.corr_util import logger
from pixcorrect import decaminfo


###########################################
def check_cti(image, CTI, verbose=0):
    """Function to check for presence of CTI"""

#
#   Initialize ctiDict
#
    ctiDict = {'isCTI': False}
    ctiDict['expnum'] = image['EXPNUM']

    # Also create the BAND and NITE keywords if they are not present
    try:
        image['BAND']
    except:
        image['BAND'] = decaminfo.get_band(image['FILTER'])
    try:
        image['NITE']
    except:
        image['NITE'] = decaminfo.get_nite(image['DATE-OBS'])

    band = image['BAND'].strip()
    sec = section2slice(image['DATASEC' + CTI['amp']])
#
#   This could become useful if it is necessary to start examining the opposite amplifier in
#   conjunction with the amplifier that is having a problem
#
#    if (CTI['amp']=="A"):
#        osec = section2slice(image['DATASEC'+'B'])
#    else:
#        osec = section2slice(image['DATASEC'+'A'])

    maxiter = 10
    converge_num = 0.0001
    clipsig = 3.0

    clip_avg, clip_med, clip_std = lb.medclip(image.data[sec], clipsig, maxiter, converge_num, verbose=0)
    logger.info(' CTI: Global(clipped): median = {:.3f}, stddev = {:.3f} '.format(clip_med, clip_std))
    ctiDict['cmed'] = float(clip_med)
    ctiDict['cstd'] = float(clip_std)
    clow = clip_med - (3.0 * clip_std)
    ctiDict['clow'] = float(clow)

#    oclip_avg,oclip_med,oclip_std=medclip(image.data[osec],clipsig,maxiter,converge_num,verbose)
#    print(" Global(oclipped): median = {:.3f}, stddev = {:.3f} ".format(oclip_med,oclip_std))
#    oclow=oclip_med-(3.0*oclip_std)

#
#   Obtain row-by-row median to look for horizontal striping (also needed to check/reject edgebleeds)
#
    row_med = np.median(image.data[sec], axis=1)
    wsm = np.where(row_med < clow)
    nrow_low = row_med[wsm].size
#
#   Hacky attempt to check for edge-bleed
#
    iedge = [4, 4091]
    while row_med[iedge[0]] < clow:
        iedge[0] = iedge[0] + 1
    while row_med[iedge[1]] < clow:
        iedge[1] = iedge[1] - 1
    if iedge[0] == 4:
        iedge[0] = 0
    if iedge[1] == 4091:
        iedge[1] = 4095
    nrow_edge = 4096 - (iedge[1] - iedge[0] + 1)
    logger.info(' CTI: Number of low rows: {:d} (nrow_edge={:d}) '.format(nrow_low, nrow_edge))

#
#   Blank out pixels that are below the 3-sigma level with respect to median
#   This removes power from vertical stripes
#
    wsm = np.where(image.data[sec] < clow)
    npix_low = image.data[sec][wsm].size
    logger.info(' CTI: Number of low pixels: {:d} '.format(npix_low))
    u = image.data[sec] - clip_med
    u[wsm] = 0.0
#
#   Harder cut currently not needed.  If used this would get rid of all pixels below the median level
#   (effectively this reduces the amount that noise suppresses contrast of the auto-correlation signal from CTI)
#
#    wsm=np.where(u<0.)
#    npix_zero=u[wsm].size
#    logger.info(' CTI: Number of sub-zero pixels: {:d} '.format(npix_zero))
#    u[wsm]=0.0

#
#   Calculate a set of auto-correlations by sampling lags in the x-direction and
#       then two diaganol sets at PA=+/-45 degrees
#   Note: y-direction lags would be succeptible to both bad columns and bleeds.
#   These are normalized by the auto-correlation with lag 0 (defined as 'a' below).
#   Take a maximum lag that will be calculated and use that to trim the image.
#   Note: This both gets rid of most edge-effects automatically but also removes the need to calculate an effective normalization for higher lags
#
    maxlag = 100
    lagList = [0, 1, 3, 5, 7, 11, 15, 19, 23, 31, 37, 45]

    a = np.sum(u[maxlag:-maxlag, maxlag:-maxlag] * u[maxlag:-maxlag, maxlag:-maxlag])
#    b=np.sum(v[maxlag:-maxlag,maxlag:-maxlag]*v[maxlag:-maxlag,maxlag:-maxlag])
    x = [1.0]
    d1 = [1.0]
    d2 = [1.0]
#    vx=[1.0]
#    vd1=[1.0]
#    vd2=[1.0]
#
#   More lags than those sampled are needed because the diagonal (PA=+/-45) measures will need to be interpolated
#   for comaparison to lags in the x-direction.
#

    for lag in lagList:
        if lag != 0:
            x.append(np.sum(u[maxlag:-maxlag, maxlag:-maxlag] * u[maxlag:-maxlag, maxlag - lag:-maxlag - lag]) / a)
            d1.append(np.sum(u[maxlag:-maxlag, maxlag:-maxlag] * u[maxlag-lag:-maxlag - lag, maxlag - lag:-maxlag - lag]) / a)
            d2.append(np.sum(u[maxlag:-maxlag, maxlag:-maxlag] * u[maxlag-lag:-maxlag - lag, maxlag + lag:-maxlag + lag]) / a)
#            vx.append(np.sum(v[maxlag:-maxlag,maxlag:-maxlag]*v[maxlag:-maxlag,maxlag-lag:-maxlag-lag])/b)
#            vd1.append(np.sum(v[maxlag:-maxlag,maxlag:-maxlag]*v[maxlag-lag:-maxlag-lag,maxlag-lag:-maxlag-lag])/b)
#            vd2.append(np.sum(v[maxlag:-maxlag,maxlag:-maxlag]*v[maxlag-lag:-maxlag-lag,maxlag+lag:-maxlag+lag])/b)

    data = {'lag': np.array(lagList),
            'x': np.array(x),
            'd1': np.array(d1),
            'd2': np.array(d2)
#          'vx':np.array(vx),
#          'vd1':np.array(vd1),
#          'vd2':np.array(vd2)
           }

    r2 = np.sqrt(2.0)
    l1 = data['lag']
    l2 = data['lag'] * r2
    x1 = data['x']
    d1i = np.interp(data['lag'], l2, data['d1'])
    d2i = np.interp(data['lag'], l2, data['d2'])
    rd1 = data['x'] / d1i
    rd2 = data['x'] / d2i

#    vx1=data['vx']
#    vd1i=np.interp(data['lag'],l2,data['vd1'])
#    vd2i=np.interp(data['lag'],l2,data['vd2'])
#    vrd1=data['vx']/vd1i
#    vrd2=data['vx']/vd2i
##    vdx=data['x']/data['vx']
#    vdx=(rd1+rd2)/(vrd1+vrd2)

    logger.info(' CTI: lags {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(l1[3], l1[4], l1[6], l1[8], l1[10]))
    logger.info(' CTI:  lx  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(x1[3], x1[4], x1[6], x1[8], x1[10]))
    logger.info(' CTI: d1i  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(d1i[3], d1i[4], d1i[6], d1i[8], d1i[10]))
    logger.info(' CTI: d2i  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(d2i[3], d2i[4], d2i[6], d2i[8], d2i[10]))
    logger.info(' CTI: ld1  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(rd1[3], rd1[4], rd1[6], rd1[8], rd1[10]))
    logger.info(' CTI: ld2  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(rd2[3], rd2[4], rd2[6], rd2[8], rd2[10]))
#    logger.info(' CTI: lvx  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vx1[3],vx1[4],vx1[6],vx1[8],vx1[10]))
#    logger.info(' CTI:vd1i  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vd1i[3],vd1i[4],vd1i[6],vd1i[8],vd1i[10]))
#    logger.info(' CTI:vd2i  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vd2i[3],vd2i[4],vd2i[6],vd2i[8],vd2i[10]))
#    logger.info(' CTI:vld1  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vrd1[3],vrd1[4],vrd1[6],vrd1[8],vrd1[10]))
#    logger.info(' CTI:vld2  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vrd2[3],vrd2[4],vrd2[6],vrd2[8],vrd2[10]))
#    logger.info(' CTI:vdx0  {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '.format(vdx[3],vdx[4],vdx[6],vdx[8],vdx[10]))

#
#   Set band dependent thresholds...
#   Note the criteria used are based on an empirical study of the one example we currently have (CCD=41, Y6)
#
    nrow_lim = 5
    if band != "Y":
        cclim = 0.9
    else:
        cclim = 1.15
#
#   Now check and set flag based on empirical critera.
#       First are the horizontal streaks that can appear...
#       Second are the comparison of the auto-correlation in the x and average of the diaganol directrions
#
    flag_cti = False
    if nrow_low - nrow_edge >= nrow_lim:
        flag_cti = True

    avg_rd = (rd1 + rd2) / 2.0
    if avg_rd[3] > cclim and avg_rd[4] > cclim and avg_rd[6] > cclim:
        flag_cti = True

    if flag_cti:
        ctiDict['isCTI'] = True

    return ctiDict


###########################################
def mask_cti(image, CTI, ctiDict, verbose=0):
    """Function to mask the region affected by a lightbulb"""

    if ctiDict['isCTI']:
        sec = section2slice(image['DATASEC' + CTI['amp']])
        image.mask[sec] |= BADPIX_BADAMP

    logger.info(' CTI: mask applied to image')

    return image
