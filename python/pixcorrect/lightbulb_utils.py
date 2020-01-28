# $Id: lightbulb_utils.py 46990 2018-05-10 19:57:25Z rgruendl $
# $Rev:: 46990                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2018-05-10 14:57:25 #$:  # Date of last commit.

"""Lightbulb masking functions
"""

import re
import numpy as np
from scipy.optimize import curve_fit
from pixcorrect.corr_util import logger
from despyfits.maskbits import *

###########################################
def check_lightbulb_explist(expnum, explist):
    """Check whether an exposure number falls within a list
        :Parameters:
            - expnum: expnum to check
            - explist: string that converts to a set of exposure numbers
    """

    if explist == "All":
        retval = True
    else:
        tmp_exp_list = explist.split(",")
        ExpList = []
        for exp in tmp_exp_list:
            if re.search("-", exp) is None:
                ExpList.append([int(exp), int(exp)])
            else:
                ExpRange = exp.split("-")
                if ExpRange[0] == '':
                    ExpList.append([1, int(ExpRange[1])])
                elif ExpRange[1] == "":
                    ExpList.append([int(ExpRange[0]), 1000000000])
                else:
                    ExpList.append([int(ExpRange[0]), int(ExpRange[1])])
        retval = False
        for ExpRange in ExpList:
            if ExpRange[0] <= expnum <= ExpRange[1]:
                retval = True

    return retval


###########################################
def medclip(data, clipsig=3.0, maxiter=10, converge_num=0.001, verbose=0):
    """ Function to examine data and determine average, median, stddev
    using a clipping algorithm
    Inputs: data: image array
            clipsig:  The number of N-sigma to be excluded when clipping
            maxiter:  Maximum number of iterations to perform
            converge_num:  Convergence Criteria

    Outputs: avgval: average from clipped distribution
             medval: median from clipped distribution
             stdval: stddev from clipped distribution
    """

    ct = data.size
    iter = 0
    c1 = 1.0
    c2 = 0.0

    avgval = np.mean(data)
    medval = np.median(data)
    sig = np.std(data)
    wsm = np.where(abs(data-medval) < clipsig * sig)
    if 0 < verbose < 4:
        logger.debug("iter,avgval,medval,sig")
    if 2 < verbose < 4:
        logger.debug("{:d} {:.2f} {:.2f} {:.2f} ".format(0, avgval, medval, sig))
    if verbose > 3:
        logger.debug("iter,avgval,medval,sig")
        logger.debug("{:d} {:.2f} {:.2f} {:.2f} {:d} {:d} {:.1f} ".format(0, avgval, medval, sig, ct, c1, c2))

    while c1 >= c2 and iter < maxiter:
        iter += 1
        lastct = ct
        avgval = np.mean(data[wsm])
        medval = np.median(data[wsm])
        sig = np.std(data[wsm])
        wsm = np.where(abs(data-medval) < clipsig * sig)
        ct = len(wsm[0])
        if ct > 0:
            c1 = abs(ct - lastct)
            c2 = converge_num * lastct
        if 2 < verbose < 4:
            logger.debug("{:d} {:.2f} {:.2f} {:.2f} ".format(iter, avgval, medval, sig))
        if verbose > 3:
            logger.debug("{:d} {:.2f} {:.2f} {:.2f} {:d} {:d} {:.1f} ".format(iter, avgval, medval, sig, ct, c1, c2))
#   End of while loop
    if iter >= maxiter:
        logger.info("Warning: medclip had not yet converged after {:d} iterations".format(iter))

    medval = np.median(data[wsm])
    avgval = np.mean(data[wsm])
    stdval = np.std(data[wsm])
    if verbose > 0:
        logger.info("{:d} {:.2f} {:.2f} {:.2f} ".format(iter + 1, avgval, medval, stdval))

    return avgval, medval, stdval


#########################################
def rgauss(r, a, b, c):
    return a * np.exp(-(r * r) / (2. * b * b)) + c

###########################################
def check_lightbulb(image, LBD, verbose=0):
    """Function to check for presence of light bulb"""

#
#   Get image statistics
#
    bulbDict = {}
    _, clip_med, clip_std = medclip(image.data, clipsig=3.0, maxiter=20, converge_num=0.001, verbose=0)
    logger.info(" LIGHTBULB: Global(clipped): median = {:.3f}, stddev = {:.3f} ".format(clip_med, clip_std))
    bulbDict['cmed'] = float(clip_med)
    bulbDict['cstd'] = float(clip_std)

#
#   Looks at central location of bulb and get some basic measure of brightness
#
    x0 = LBD['xc']
    y0 = LBD['yc']
    iy1 = y0 - 10
    iy2 = y0 + 10
    ix1 = x0 - 10
    ix2 = x0 + 10
    central_med = np.median(image.data[iy1:iy2, ix1:ix2])
    bulbDict['bulb_sb'] = float(central_med) - bulbDict['cmed']
    logger.info(" LIGHTBULB: Potential Bulb Central SB: median = {:.3f}".format(bulbDict['bulb_sb']))

    BulbSig = (central_med - clip_med) / clip_std
    bulbDict['bulb_sig'] = (float(BulbSig) * 21.)
    logger.info(" LIGHTBULB: Bulb Significance = {:.3f} sigma.".format(BulbSig * 21.))
    central_sat = np.sum(image.mask[iy1:iy2, ix1:ix2] & 2)/ 2
    logger.info(" LIGHTBULB: Number of central saturated pixels: {:d} ".format(int(central_sat)))
    bulbDict['num_sat'] = int(central_sat)
#
#   Form a radial profile and then fit with Gaussian centered at r=0
#
#    image.data[np.isnan(image.data)]=0
    y, x = np.indices(image.data.shape)
    r = np.sqrt((x - x0)**2 + (y - y0)**2)

    rbinsize = 10.0
    radbin = np.arange(0., LBD['rad'], rbinsize)
    radbin_c = (radbin[0:-1] + radbin[1:]) / 2.0
    medbin = np.zeros(radbin_c.size)
    stdbin = np.zeros(radbin_c.size)
    for i in range(radbin_c.size):
        wsm = np.where(np.logical_and(r >= radbin[i], r < radbin[i + 1]))
        if image.data[wsm].size > 0:
            medbin[i] = np.median(image.data[wsm])
            stdbin[i] = np.std(image.data[wsm])
            if medbin[i] > bulbDict['cmed']:
                stdbin[i] = np.sqrt(medbin[i] - bulbDict['cmed'] + 36.0)
            else:
                stdbin[i] = np.sqrt(36.0)
        else:
            medbin[i] = 0.0
            stdbin[i] = np.sqrt(36.0)
#
#   Make a simple estimate of radius to mask (then confirm with Gaussian fit)
#
    bkg_limit = bulbDict['cmed'] + bulbDict['cstd']
    i = 0
    while medbin[i] > bkg_limit and i < radbin_c.size - 1:
        i += 1
    logger.info(" LIGHTBULB: Median Radial Profile estimate shows level of {:.2f} at a radius of {:.2f} ".format(medbin[i], radbin_c[i]))
    bulbDict['rad'] = float(radbin_c[i])
#
#   Guess for FWHM, Amp, bkg
#
    i = 0
    while medbin[i] - bulbDict['cmed'] > 0.5 * bulbDict['bulb_sb'] and i < radbin_c.size - 1:
        i += 1
    sig_guess = radbin_c[i] / 2.355
    amp_guess = bulbDict['bulb_sb'] / (sig_guess * 2.355)
    r_guess = [amp_guess, sig_guess, bulbDict['cmed']]
#
#   Fit radial profile
#
    try:
        popt, pcov = curve_fit(rgauss, radbin_c, medbin, p0=r_guess, sigma=stdbin, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        logger.info(" LIGHTBULB: Gauss FIT results (amp, sig, bkg): {:.2f} {:.2f} {:.2f} ".format(popt[0], popt[1], popt[2]))
        logger.info(" LIGHTBULB: Gauss FIT results     perr(covar): {:.2f} {:.2f} {:.2f} ".format(perr[0], perr[1], perr[2]))
        #g_fitval = rgauss(radbin_c, *popt)
        bulbDict['g_amp'] = popt[0]
        bulbDict['g_wid'] = popt[1] * 2.355
        bulbDict['g_bkg'] = popt[2]
        bulbDict['g_amperr'] = perr[0]
        bulbDict['g_widerr'] = perr[1] * 2.355
        bulbDict['g_bkgerr'] = perr[2]
    except:
        logger.info(" LIGHTBULB: Gauss FIT FAILED")
#        popt=(0.0,1.0,bulbDict['cmed'])
        bulbDict['g_amp'] = -1.0
        bulbDict['g_wid'] = -1.0
        bulbDict['g_bkg'] = -1.0
        bulbDict['g_amperr'] = -1.0
        bulbDict['g_widerr'] = -1.0
        bulbDict['g_bkgerr'] = -1.0

    return bulbDict


###########################################
def mask_lightbulb(image, LBD, bulbDict, verbose=0):
    """Function to mask the region affected by a lightbulb"""

#
#   Empirically derived relation for the width of the columnar masks
#       pw1 --> width of region extending toward the read registers
#       pw2 --> width of region extending away from the read registers
#
    if bulbDict['bulb_sb'] > 0.:
        y1 = np.power(bulbDict['bulb_sb'], 0.75) / 10.
        y2 = np.power(bulbDict['bulb_sb'], 0.75) / 100.
        pw1 = 150. * (np.log10(y1) - 0.5)
        pw2 = 150. * (np.log10(y2) - 0.5)
        if pw1 < 16.0:
            pw1 = 0.0
        if pw2 < 16.0:
            pw2 = 0.0
    else:
        pw1 = 0.0
        pw2 = 0.0

    xb0 = LBD['xc']
    yb0 = LBD['yc']
    yb, xb = np.indices(image.mask.shape)
    rb = np.sqrt((xb - xb0)**2 + (yb - yb0)**2)

    rlim = 1.66 * bulbDict['rad']
    xb11 = xb0 - pw1
    xb12 = xb0 + pw1
    xb21 = xb0 - pw2
    xb22 = xb0 + pw2

    if pw2 > 7.:
        wsm3 = np.where(np.logical_or(rb < rlim, np.logical_or(np.logical_and(yb >= yb0, np.logical_and(xb > xb11, xb < xb12)),
                                                               np.logical_and(yb <= yb0, np.logical_and(xb > xb21, xb < xb22)))))
        logger.info(" LIGHTBULB: masking circular region with radius {:.1f} ".format(rlim))
        logger.info(" LIGHTBULB: masking columnar region for y>{:.1f} with width {:.1f} ".format(yb0, pw1 * 2.))
        logger.info(" LIGHTBULB: masking columnar region for y<{:.1f} with width {:.1f} ".format(yb0, pw2 * 2.))
    elif pw1 > 7.:
        wsm3 = np.where(np.logical_or(rb < rlim, np.logical_and(yb >= yb0, np.logical_and(xb > xb11, xb < xb12))))
        logger.info(" LIGHTBULB: masking circular region with radius {:.1f} ".format(rlim))
        logger.info(" LIGHTBULB: masking columnar region for y>{:.1f} with width {:.1f} ".format(yb0, pw1 * 2.))
    else:
        wsm3 = np.where(rb < rlim)
        logger.info(" LIGHTBULB: masking circular region with radius {:.1f} ".format(rlim))

    image.mask[wsm3] |= BADPIX_BPM

    logger.info(" LIGHTBULB: mask applied to image")

    return image
