#!/usr/bin/env python3

# $Id: connect_streaks.py 47132 2018-06-13 19:05:15Z rgruendl $
# $Rev:: 47132                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2018-06-13 14:05:15 #$:  # Date of last commit.
#
# Note original code was developed by Gary Bernstein...
#

"""
Read streak tables, determine candidate missed streaks, add them to masks
"""

import re
import time
import shutil
import numpy as np
from scipy.optimize import newton
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as pl

import fitsio
from despyastro import wcsutil
from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage
from pixcorrect.PixCorrectDriver import PixCorrectDriver, filelist_to_list
from pixcorrect.corr_util import logger

# Which section of the config file to read for this step
config_section = 'connect_streaks'

# Utility routines for this task:
def gnomonic(ra, dec, ra0, dec0):
    # ra, dec are the target points
    # ra0, dec0 are the projection centers
    # All are given in degrees
    # returns x, y also in degrees, with
    # x to east and y to north
    dtor = np.pi / 180.
    cd = np.cos(dec * dtor)
    cd0 = np.cos(dec0 * dtor)
    sd = np.sin(dec * dtor)
    sd0 = np.sin(dec0 * dtor)
    cr = np.cos((ra - ra0) * dtor)
    sr = np.sin((ra - ra0) * dtor)
    x = (cd * sr) / (sd0 * sd + cd0 * cd * cr) / dtor
    y = (cd0 * sd - sd0 * cd * cr) / (sd0 * sd + cd0 * cd * cr) / dtor
    return x, y

# And its inverse
def gnomonicInverse(x, y, ra0, dec0):
    # take location (x,y) in
    # gnomonic projection about ra0,dec0
    # and return ra,dec.
    # Everything in degrees.
    dtor = np.pi / 180.
    rho = np.hypot(x, y) * dtor
    c = np.arctan(rho)
    cd0 = np.cos(dec0 * dtor)
    sd0 = np.sin(dec0 * dtor)
    dec = np.arcsin(np.cos(c) * sd0 + y * dtor * np.sin(c) * cd0 / rho) / dtor
    ra = np.arctan(x * dtor * np.sin(c) / (rho * cd0 * np.cos(c) - y * dtor * sd0 * np.sin(c))) / dtor + ra0
    return ra, dec

def boxCross(x0, y0, mx_in, my_in, xmin, xmax, ymin, ymax):
    # Function that determines what (if any) parts of the
    # lines defined by (x-x0)*my = (y-y0)*mx cross
    # the rectangle bounded by xmin<x<xmax and
    # ymin<y<ymax.  Returns four arrays,
    # x1,x2,y1,y2 giving endpoints of crossings.
    # x2<x1 signifies no intersection.

    # First standardize lines to have mx>0.  If mx==0,
    # special handling needed.
    ss = np.sign(mx_in)
    mx = mx_in * ss
    my = np.where(ss < 0, -my_in, my_in)

    # Find the two points where line crosses x=xmin,x=xmax
    xcross1 = np.where(ss == 0, x0, xmin)
    ycross1 = np.where(ss == 0, ymin, y0 + my * (xmin - x0) / mx)
    xcross2 = np.where(ss == 0, x0, xmax)
    ycross2 = np.where(ss == 0, ymax, y0 + my * (xmax - x0) / mx)

    # No intersection if y does not enter range at all
    miss = np.logical_or(np.logical_and(ycross1 < ymin, ycross2 < ymin),
                         np.logical_and(ycross1 > ymax, ycross2 > ymax))

    # or if it's a vertical line and out of x range
    miss = np.logical_or(miss, np.logical_and(ss == 0, x0 < xmin))
    miss = np.logical_or(miss, np.logical_and(ss == 0, x0 > xmax))
    # or if it's a horizontal line and out of y range
    miss = np.logical_or(miss, np.logical_and(my == 0, y0 < ymin))
    miss = np.logical_or(miss, np.logical_and(my == 0, y0 > ymax))

    # Pull crossing points to be within [ymin,ymax]
    fix = ycross1 > ymax
    xcross1 = np.where(fix, x0 + (ymax - y0) * mx / my, xcross1)
    ycross1 = np.where(fix, ymax, ycross1)
    fix = ycross2 > ymax
    xcross2 = np.where(fix, x0 + (ymax - y0) * mx / my, xcross2)
    ycross2 = np.where(fix, ymax, ycross2)
    fix = ycross1 < ymin
    xcross1 = np.where(fix, x0 + (ymin - y0) * mx / my, xcross1)
    ycross1 = np.where(fix, ymin, ycross1)
    fix = ycross2 < ymin
    xcross2 = np.where(fix, x0 + (ymin - y0) * mx / my, xcross2)
    ycross2 = np.where(fix, ymin, ycross2)

    # Backfill non-crossers
    xcross1[miss] = 0
    xcross2[miss] = -1.
    ycross1[miss] = 0
    ycross2[miss] = 0

    # Special case for valid horizontal lines?? (my==0)
    fix = np.logical_and(my == 0, y0 > ymin)
    fix = np.logical_and(fix, y0 < ymax)
    xcross1 = np.where(fix, xmin, xcross1)
    xcross2 = np.where(fix, xmax, xcross2)
    ycross1 = np.where(fix, y0, ycross1)
    ycross2 = np.where(fix, y0, ycross2)
    return xcross1, xcross2, ycross1, ycross2

class Line:
    '''
    Class representing the line connecting two points.
    Used to establish a new (t,u) coordinate system
    where t is along the line and u is perpendicular distance.
    Inputs to all functions can be scalars or numpy arrays.
    '''
    def __init__(self, x1, y1, x2, y2):
        # Define line connecting (x1,y1) to (x2,y2)
        self.x0 = 0.5 * (x1 + x2)
        self.y0 = 0.5 * (y1 + y2)
        dx = x2 - x1
        dy = y2 - y1
        self.mx = dx / np.hypot(dx, dy)
        self.my = dy / np.hypot(dx, dy)

    def xy2tu(self, x, y):
        # Transform x,y coordinates to the (t,u) system defined by line
        dx = x - self.x0
        dy = y - self.y0
        return self.mx * dx + self.my * dy, -self.my * dx + self.mx * dy

    def tu2xy(self, t, u):
        # Transform t,u coordinates in the line's system into x,y.
        dx = self.mx * t - self.my * u
        dy = self.my * t + self.mx * u
        return dx + self.x0, dy + self.y0

    def ux2t(self, u, x):
        # Solve for the t value where u and x are given
        return (self.my * u + (x - self.x0)) / self.mx

    def uy2t(self, u, y):
        # Solve for the t value where u and y are given
        return (-self.mx * u + (y - self.y0)) / self.my

def boxTrack(line, w, xmin, xmax, ymin, ymax):
    ''' Determine the corners of the minimal rectangle centered
    on the line and having width 2w that contains all possible points
    rectangle bounded by [xy][min|max].
        Return value will be 4x2 array of these corners' xy coords.  They
    will all be zero if there is no intersection.
    '''
    xbox = np.array([xmin, xmax, xmax, xmin, xmin])
    ybox = np.array([ymin, ymin, ymax, ymax, ymin])
    xVaries = (True, False, True, False)  # which coord varies along edge?
    fixedValue = (ymin, xmax, ymax, xmin) # what is value of fixed x or y?
    t, u = line.xy2tu(xbox, ybox)
    hw = 0.5 * w
    # Make an array which is +1,0,-1 as the corner is
    # above, inside, or below the track of width w
    corner_state = np.where(u > hw, 1, 0)
    corner_state = np.where(u < hw, -1, corner_state)

    # Analyze each corner and edge of the box
    # to find t values of all crossings of track
    # with box edges, or of box corners in the track

    # All corners within track are possible extrema
    t_extremes = t[corner_state == 0].tolist()

    for corner in range(4):
        state1 = corner_state[corner]
        state2 = corner_state[corner + 1]
        if (state1 > 0 and state2 <= 0) or (state2 > 0 and state1 <= 0):
            # There should be an upper crossing on this segment
            if xVaries[corner]:
                t_extremes.append(line.uy2t(+hw, fixedValue[corner]))
            else:
                t_extremes.append(line.ux2t(+hw, fixedValue[corner]))
        if (state1 < 0 and state2 >= 0) or (state2 < 0 and state1 >= 0):
            # There should be a lower crossing on this segment
            if xVaries[corner]:
                t_extremes.append(line.uy2t(-hw, fixedValue[corner]))
            else:
                t_extremes.append(line.ux2t(-hw, fixedValue[corner]))
    # Now the limits of t for the rectangle will be min and max of box
    out = np.zeros((4, 2), dtype=float)
    if t_extremes:
        # Only do this if there is any overlap
        tmin = np.min(t_extremes)
        tmax = np.max(t_extremes)
        rectangle_t = np.array([tmin, tmin, tmax, tmax])
        rectangle_u = np.array([-hw, +hw, +hw, -hw]) # ?? backwards ??
        rectangle_x, rectangle_y = line.tu2xy(rectangle_t, rectangle_u)
        out[:, 0] = rectangle_x
        out[:, 1] = rectangle_y
    return out

def friends(xc1, xc2, yc1, yc2, mx, my, ccdnum, i, j, max_sine=0.02):
    # Determine whether the two streaks
    # pass within tolerance of each other's centers

    # Same-CCD streaks cannot be friends
    # (Though they may get linked if both friends of
    # a streak on another CCD.)
    if ccdnum[i] == ccdnum[j]:
        return False

    xi = (xc1[i] + xc2[i]) * 0.5
    xj = (xc1[j] + xc2[j]) * 0.5
    yi = (yc1[i] + yc2[i]) * 0.5
    yj = (yc1[j] + yc2[j]) * 0.5

    di = (xi - xj) * my[j] - (yi - yj) * mx[j]
    dj = (xj - xi) * my[i] - (yj - yi) * mx[i]

    # This criterion is that the sine of angle between the
    # center-to-center line and each streak is < max_sine
    dij = np.hypot(xi - xj, yi - yj)
    return max(np.abs(di), np.abs(dj)) < max_sine * dij


class ConnectStreaks(PixCorrectDriver):
    description = "Predict missed streak detections and mask them"
    step_name = config_section

    DEFAULT_STREAK_NAME_IN = 'streak'
    DEFAULT_STREAK_NAME_OUT = 'streak2'
    DEFAULT_IMAGE_NAME_IN = 'hmmasked'
    DEFAULT_IMAGE_NAME_OUT = 'immasked'
    DEFAULT_ADD_WIDTH = 0
    DEFAULT_MAX_EXTRAPOLATE = 1.1 * np.hypot(2048., 4096.) * 0.263 / 3600 # A bit more than 1 CCD diagonal

    @classmethod
    def __call__(cls, streak_list, image_list,
                 streak_name_in, streak_name_out,
                 image_name_in, image_name_out,
                 add_width, max_extrapolate,
                 plotfile=None):

        """
        Read input list of streak detections and predict where a streak
        crossed a CCD but was missed.  Then create new copies of images,
        altering masks to set STREAK bit in new streaks.

        :Parameters:
            - `streak_list`: list of input streak file names
            - `image_list`: list of names of image files to be updated
            - `streak_name_in`: string to replace in input streak filenames
            - `streak_name_out`: replacement string for output streak filenames
            - `image_name_in`: string to replace in input image filenames
            - `image_name_out`: replacement string for output image filenames
            - `add_width`:  number of pixels to grow (or shrink) streak width
            - `max_extrapolate`: farthest to start a new streak from endpoint of an existing one (degrees)
            - `plotfile`: if given, a diagram of streaks is drawn into this file
        """

        logger.info('Reading {:d} streak files'.format(len(streak_list)))

        # Read in all the streak RA/Dec, into a dictionary keyed by CCDNUM,
        # which should be in the primary header.  Also save a dictionary of
        # the file names for these
        streak_corners = {}
        streak_names = {}
        for streakfile in streak_list:
            logger.info(f"Reading streak file {streakfile}")
            with fitsio.FITS(streakfile, 'r') as fits:
                ccdnum = fits[0].read_header()['CCDNUM']
                streak_names[ccdnum] = streakfile
                tab = fits[1].read()
                if len(tab) > 0:
                    streak_corners[ccdnum] = fits[1].read()['CORNERS_WCS']

        logger.info('Reading WCS from {:d} CCDs'.format(len(image_list)))

        # Read in the WCS for each CCD for which we have an image,
        # also put into dicts keyed by CCDNUM
        # Will get these directly from FITS instead of using DESImage in order
        # to save reading all of the data.
        wcs = {}
        crval1 = []
        crval2 = []
        for imgfile in image_list:
            try:
                hdr = fitsio.read_header(imgfile, 0)
                ccd = hdr['CCDNUM']
                crval1.append(hdr['CRVAL1'])
                crval2.append(hdr['CRVAL2'])
                # Due to a bug in fitsio 1.0.0rc1+0, we need to clean up the
                # header before feeding it to wcsutil and remove the 'None' and other problematic items
                for k in hdr:
                    # Try to access the item, if failed we have to remove it
                    if not k:
                        hdr.delete(k)
                        continue
                    try:
                        _ = hdr[k]
                    except:
                        logger.info("Removing keyword: {:s} from header".format(k))
                        hdr.delete(k)
                wcs[ccd] = wcsutil.WCS(hdr)
            except Exception as e:
                print(e) ###
                logger.error('Failure reading WCS from {:s}'.format(imgfile))
                return 1

        # Determine a center for local gnomonic projection
        ra0 = np.median(crval1)
        dec0 = np.median(crval2)

        # Calculate upper and lower bounds of each CCD in the local
        # gnomonic system.
        ccd_x1 = np.zeros(63, dtype=float)
        ccd_x2 = np.zeros(63, dtype=float)
        ccd_y1 = np.zeros(63, dtype=float)
        ccd_y2 = np.zeros(63, dtype=float)

        ccd_xmin = 1.
        ccd_xmax = 2048.
        ccd_ymin = 1.
        ccd_ymax = 4096.
        ccd_corners_xpix = np.array([ccd_xmin, ccd_xmin, ccd_xmax, ccd_xmax])
        ccd_corners_ypix = np.array([ccd_ymin, ccd_ymax, ccd_ymax, ccd_ymin])
        for ccd, w in wcs.items():
            ra, dec = w.image2sky(ccd_corners_xpix, ccd_corners_ypix)
            x_corners, y_corners = gnomonic(ra, dec, ra0, dec0)
            ccd_x1[ccd] = np.min(x_corners)
            ccd_y1[ccd] = np.min(y_corners)
            ccd_x2[ccd] = np.max(x_corners)
            ccd_y2[ccd] = np.max(y_corners)

        # Now collect information on all of the streak segments that we have
        ccdnum = []
        ra_corner = []
        dec_corner = []

        for ccd, streaks in streak_corners.items():
            if ccd not in wcs:
                # Skip segments on CCDs that have no WCS
                logger.warning('No WCS found for streaks on CCD {:d}'.format(ccd))
                continue
            n1, _, _ = streaks.shape
            for i in range(n1):
                ccdnum.append(ccd)
                ra_corner.append(streaks[i, :, 0])
                dec_corner.append(streaks[i, :, 1])
        # Put streak corners into gnomonic system for this exposure
        x1, y1 = gnomonic(np.array([r[0] for r in ra_corner], dtype=float),
                          np.array([d[0] for d in dec_corner], dtype=float),
                          ra0, dec0)
        x2, y2 = gnomonic(np.array([r[1] for r in ra_corner], dtype=float),
                          np.array([d[1] for d in dec_corner], dtype=float),
                          ra0, dec0)
        x3, y3 = gnomonic(np.array([r[2] for r in ra_corner], dtype=float),
                          np.array([d[2] for d in dec_corner], dtype=float),
                          ra0, dec0)
        x4, y4 = gnomonic(np.array([r[3] for r in ra_corner], dtype=float),
                          np.array([d[3] for d in dec_corner], dtype=float),
                          ra0, dec0)
        ccdnum = np.array(ccdnum, dtype=int)

        # Describe each segmet by two endpoints at the midpoints of short sides
        # Will need to decide which is the short side
        d12 = np.hypot(x2 - x1, y2 - y1)
        d23 = np.hypot(x3 - x2, y3 - y2)
        xleft = np.where(d12 < d23, 0.5 * (x1 + x2), 0.5 * (x2 + x3))
        yleft = np.where(d12 < d23, 0.5 * (y1 + y2), 0.5 * (y2 + y3))
        xright = np.where(d12 < d23, 0.5 * (x3 + x4), 0.5 * (x4 + x1))
        yright = np.where(d12 < d23, 0.5 * (y3 + y4), 0.5 * (y4 + y1))
        dx = xright - xleft
        dy = yright - yleft
        # Calculate a width as 2x the
        # largest perp distance from a vertex to this line
        w1 = np.abs(dx * (y1 - yleft) - dy * (x1 - xleft)) / np.hypot(dx, dy)
        w2 = np.abs(dx * (y2 - yleft) - dy * (x2 - xleft)) / np.hypot(dx, dy)
        w3 = np.abs(dx * (y3 - yleft) - dy * (x3 - xleft)) / np.hypot(dx, dy)
        w4 = np.abs(dx * (y4 - yleft) - dy * (x4 - xleft)) / np.hypot(dx, dy)
        wmax = np.maximum(w1, w2)
        wmax = np.maximum(wmax, w3)
        wmax = np.maximum(wmax, w4)
        wmax = 2 * wmax

        # Rearrange so that xleft <= xright
        swapit = xright < xleft
        tmp = np.where(swapit, xleft, xright)
        xleft = np.where(swapit, xright, xleft)
        xright = np.array(tmp)
        tmp = np.where(swapit, yleft, yright)
        yleft = np.where(swapit, yright, yleft)
        yright = np.array(tmp)

        # Get the crossing points of the lines into CCDs
        xc1, xc2, yc1, yc2 = boxCross(xleft, yleft, dx, dy,
                                      ccd_x1[ccdnum], ccd_x2[ccdnum], ccd_y1[ccdnum], ccd_y2[ccdnum])

        # Get rid of segments that appear to miss their host CCDs
        miss = xc2 < xc1

        # Take 1st crossing point instead of left point if it has higher x, or vertical
        # with higher y, i.e. truncate the track segment at the edge of the CCD.
        replace = np.where(dx == 0, yc1 > yleft, xc1 > xleft)
        xc1 = np.where(replace, xc1, xleft)
        yc1 = np.where(replace, yc1, yleft)
        # Likewise truncate segment at right-hand crossing
        replace = np.where(dx == 0, yc2 < yright, xc2 < xright)
        xc2 = np.where(replace, xc2, xright)
        yc2 = np.where(replace, yc2, yright)

        # Backfill the non-intersections again - note that above
        # maneuvers will leave xc2<xc1 for streaks that miss their CCDs,
        # unless vertical ???
        xc1[miss] = 0.
        xc2[miss] = -1.

        # Get a final verdict on hit or miss
        miss = np.where(dx == 0, yc2 < yc1, xc2 < xc1)

        # Save information on all valid streaks
        xc1 = xc1[~miss]
        xc2 = xc2[~miss]
        yc1 = yc1[~miss]
        yc2 = yc2[~miss]
        wmax = wmax[~miss]
        ccdnum = ccdnum[~miss]

        # Express segments as slopes and midpoints
        dx = xc2 - xc1
        dy = yc2 - yc1
        mx = dx / np.hypot(dx, dy)
        my = dy / np.hypot(dx, dy)

        # Mark segments that are probably spurious edge detections
        EDGE_SLOPE = 2.  # Degrees from horizontal for edge streaks
        EDGE_DISTANCE = 0.005 # Max degrees from streak center to CCD edge for spurious streaks
        horizontal = np.abs(my) < np.sin(EDGE_SLOPE * np.pi / 180.)
        ymid = 0.5 * (yc1 + yc2)
        nearedge = np.logical_or(ccd_y2[ccdnum] - ymid < EDGE_DISTANCE,
                                 ymid-ccd_y1[ccdnum] < EDGE_DISTANCE)
        nearedge = np.logical_and(nearedge, horizontal)

        # Check short edges too
        vertical = np.abs(mx) < np.sin(EDGE_SLOPE * np.pi / 180.)
        xmid = 0.5 * (xc1 + xc2)
        tmp = np.logical_or(ccd_x2[ccdnum] - xmid < EDGE_DISTANCE,
                            xmid-ccd_x1[ccdnum] < EDGE_DISTANCE)
        nearedge = np.logical_or(nearedge, np.logical_and(tmp, vertical))

        # Decide which segments are "friends" of each other.
        # To be a friend, the center of each must be close
        # to the extension of the line of the other.
        # Accumulate a list of tracks, each track is a list of
        # individual streaks that are friends of friends
        tracks = []

        for i in range(len(xc1)):
            if nearedge[i]:
                continue  # Do not use edge tracks
            itstrack = [i] # start new track with just this
            for t in tracks:
                # Search other tracks for friends
                for j in t:
                    if friends(xc1, xc2, yc1, yc2, mx, my, ccdnum, i, j):
                        itstrack += t   # Merge track
                        tracks.remove(t) # Get rid of old one
                        break           # No need to check others
            tracks.append(itstrack)

        # Now iterate through tracks, seeing if they have missing segments
        # Create arrays to hold information on new tracks
        new_ccdnum = []
        new_xc1 = []
        new_xc2 = []
        new_yc1 = []
        new_yc2 = []
        new_ra1 = []
        new_ra2 = []
        new_dec1 = []
        new_dec2 = []
        new_width = []
        new_extrapolated = []
        new_nearest = []

        for t in tracks:
            if len(t) < 2:
                continue # Do not extrapolate singlet tracks
            ids = np.array(t)  # Make an array of indices of segments in this track
            # Fit a quadratic path to the streak endpoints
            xx = np.concatenate((xc1[ids], xc2[ids]))
            yy = np.concatenate((yc1[ids], yc2[ids]))

            # If the track slope is mostly along x, then we'll have the independent
            # variable xx be x and dependent yy will be y.  But if track
            # is more vertical, then we'll look at functions x(y) instead.
            xOrder = np.median(np.abs(mx[ids])) > np.median(np.abs(my[ids]))
            if not xOrder:
                xx, yy = yy, xx

            # Record limits of detected tracks' independent variable
            xxmin = np.min(xx)
            xxmax = np.max(xx)

            # Fit a quadratic to the points, or
            # linear if only one streak
            # Allow up to nclip points to clip
            RESID_TOLERANCE = 6. / 3600. # Clip >6" deviants
            nclip = 2
            for i in range(nclip+1):
                if len(xx) > 2:
                    A = np.vstack((np.ones_like(xx), xx, xx * xx))
                else:
                    A = np.vstack((np.ones_like(xx), xx))
                coeffs = np.linalg.lstsq(A.T, yy)[0]
                resid = yy - np.dot(A.T, coeffs)
                j = np.argmax(np.abs(resid))
                if i == nclip or np.abs(resid[j]) < RESID_TOLERANCE:
                    break
                xx = np.delete(xx, j)
                yy = np.delete(yy, j)

            # Calculate the y(x1),y(x2) where tracks
            # cross the left/right of every CCD, then
            # find the ones that will cross CCD's y.

            # These are CCD bounds, with xx being the quadratic's argument
            if xOrder:
                xx1 = ccd_x1
                xx2 = ccd_x2
                yy1 = ccd_y1
                yy2 = ccd_y2
            else:
                xx1 = ccd_y1
                xx2 = ccd_y2
                yy1 = ccd_x1
                yy2 = ccd_x2

            if len(coeffs) == 2:
                A2 = np.vstack((np.ones_like(xx2), xx2)).T
                A1 = np.vstack((np.ones_like(xx1), xx1)).T
            else:
                A2 = np.vstack((np.ones_like(xx2), xx2, xx2 * xx2)).T
                A1 = np.vstack((np.ones_like(xx1), xx1, xx1 * xx1)).T

            # yyc[12] are the dependent coordinate at crossings of xx[12] bounds
            yyc1 = np.dot(A1, coeffs)
            yyc2 = np.dot(A2, coeffs)
            # Now we ask whether the y value of streak at either edge crossing
            # is in the y range of a CCD
            missed = np.logical_or(np.maximum(yyc1, yyc2) < yy1, np.minimum(yyc1, yyc2) > yy2)
            # Also skip any CCD where we already have a streak
            for iccd in ccdnum[ids]:
                missed[iccd] = True
            missed[0] = True  # There is no CCD0
            missed[61] = True # Never use this one either, it's always dead

            # Now find intersection of new streaks with edges of their CCDs
            # Define a function for the streak path that we'll use for solving
            def poly(x, coeffs, ysolve):
                y = coeffs[0] + x * coeffs[1]
                if len(coeffs) > 2:
                    y += coeffs[2] * x * x
                return y - ysolve

            EDGE_TOLERANCE = 0.2 / 3600.  # Find x/y of edge to this accuracy (0.2 arcsec)
            for iccd in np.where(~missed)[0]:
                # This is a loop over every CCD that the track crosses but has no detected segment
                # Determine an (xx,yy) pair for its entry and exit from the CCD
                new_yy1 = yyc1[iccd]
                new_yy2 = yyc2[iccd]
                new_xx1 = xx1[iccd]
                new_xx2 = xx2[iccd]
                # left side:
                if new_yy1 < yy1[iccd]:
                    new_xx1 = newton(poly, new_xx1, args=(coeffs, yy1[iccd]), tol=EDGE_TOLERANCE)
                elif new_yy1 > yy2[iccd]:
                    new_xx1 = newton(poly, new_xx1, args=(coeffs, yy2[iccd]), tol=EDGE_TOLERANCE)
                new_yy1 = poly(new_xx1, coeffs, 0.)
                # right side
                if new_yy2 < yy1[iccd]:
                    new_xx2 = newton(poly, new_xx2, args=(coeffs, yy1[iccd]), tol=EDGE_TOLERANCE)
                elif new_yy2 > yy2[iccd]:
                    new_xx2 = newton(poly, new_xx2, args=(coeffs, yy2[iccd]), tol=EDGE_TOLERANCE)
                new_yy2 = poly(new_xx2, coeffs, 0.)
                # Does the solution lie outside the input streaks?
                extrapolated = new_xx1 < xxmin or new_xx2 > xxmax
                width = np.median(wmax[ids])

                # Calculate distance to nearest unclipped streak member
                nearest = min(np.min(np.hypot(xx - new_xx1, yy - new_yy1)),
                              np.min(np.hypot(xx - new_xx2, yy - new_yy2)))

                if not xOrder:
                    # swap xx,yy back if we had y as the independent variable
                    new_xx1, new_yy1 = new_yy1, new_xx1
                    new_xx2, new_yy2 = new_yy2, new_xx2

                # Project the coordinates back to RA, Dec
                ra1, dec1 = gnomonicInverse(new_xx1, new_yy1, ra0, dec0)
                ra2, dec2 = gnomonicInverse(new_xx2, new_yy2, ra0, dec0)

                # Append this streak to list of new ones
                new_ccdnum.append(iccd)
                new_xc1.append(new_xx1)
                new_xc2.append(new_xx2)
                new_yc1.append(new_yy1)
                new_yc2.append(new_yy2)
                new_ra1.append(ra1)
                new_ra2.append(ra2)
                new_dec1.append(dec1)
                new_dec2.append(dec2)
                new_width.append(width)
                new_extrapolated.append(extrapolated)
                new_nearest.append(nearest)

        # Make all lists into arrays
        new_ccdnum = np.array(new_ccdnum, dtype=int)
        new_xc1 = np.array(new_xc1, dtype=float)
        new_xc2 = np.array(new_xc2, dtype=float)
        new_yc1 = np.array(new_yc1, dtype=float)
        new_yc2 = np.array(new_yc2, dtype=float)
        new_ra1 = np.array(new_ra1, dtype=float)
        new_ra2 = np.array(new_ra2, dtype=float)
        new_dec1 = np.array(new_dec1, dtype=float)
        new_dec2 = np.array(new_dec2, dtype=float)
        new_width = np.array(new_width, dtype=float)
        new_extrapolated = np.array(new_extrapolated, dtype=bool)
        new_nearest = np.array(new_nearest, dtype=float)

        # Decide which new segments will be masked
        maskit = np.logical_or(~new_extrapolated, new_nearest <= max_extrapolate)

        logger.info('Identified {:d} missing streak segments for masking'.format(\
                    np.count_nonzero(maskit)))

        # Make the diagnostic plot if desired
        if plotfile is not None:
            pl.figure(figsize=(6, 6))
            pl.xlim(-1.1, 1.1)
            pl.ylim(-1.1, 1.1)
            pl.gca().set_aspect('equal')

            # Draw CCD outlines and numbers
            for ccd, w in wcs.items():
                ra, dec = w.image2sky(ccd_corners_xpix, ccd_corners_ypix)
                x_corners, y_corners = gnomonic(ra, dec, ra0, dec0)
                x = x_corners.tolist()
                y = y_corners.tolist()
                x.append(x[0])
                y.append(y[0])
                pl.plot(x, y, 'k-', label=None)
                x = np.mean(x_corners)
                y = np.mean(y_corners)
                pl.text(x, y, str(ccd), horizontalalignment='center',
                        verticalalignment='center', fontsize=14)



            # Draw input streaks marked as edge
            labelled = False
            for i in np.where(nearedge)[0]:
                x = (xc1[i], xc2[i])
                y = (yc1[i], yc2[i])
                if not labelled:
                    pl.plot(x, y, 'm-', lw=2, label='edge')
                    labelled = True
                else:
                    pl.plot(x, y, 'm-', lw=2, label=None)

            # Draw linked tracks
            s = set()
            for t in tracks:
                if len(t) > 1:
                    s = s.union(set(t))
            labelled = False
            for i in s:
                x = (xc1[i], xc2[i])
                y = (yc1[i], yc2[i])
                if not labelled:
                    pl.plot(x, y, 'b-', lw=2, label='connected')
                    labelled = True
                else:
                    pl.plot(x, y, 'b-', lw=2, label=None)

            # Draw singleton tracks as those that are neither edge nor connected
            s = s.union(set(np.where(nearedge)[0]))
            single = set(range(len(xc1)))
            single = single.difference(s)
            labelled = False
            for i in single:
                x = (xc1[i], xc2[i])
                y = (yc1[i], yc2[i])
                if not labelled:
                    pl.plot(x, y, 'c-', lw=2, label='unconnected')
                    labelled = True
                else:
                    pl.plot(x, y, 'c-', lw=2, label=None)

            # Draw missed tracks that will be masked
            labelled = False
            for i in np.where(maskit)[0]:
                x = (new_xc1[i], new_xc2[i])
                y = (new_yc1[i], new_yc2[i])
                if not labelled:
                    pl.plot(x, y, 'r-', lw=2, label='new masked')
                    labelled = True
                else:
                    pl.plot(x, y, 'r-', lw=2, label=None)


            # Draw missed tracks that will not be masked
            labelled = False
            for i in np.where(~maskit)[0]:
                x = (new_xc1[i], new_xc2[i])
                y = (new_yc1[i], new_yc2[i])
                if not labelled:
                    pl.plot(x, y, 'r:', lw=2, label='new skipped')
                    labelled = True
                else:
                    pl.plot(x, y, 'r:', lw=2, label=None)

            # legend
            pl.legend(framealpha=0.3, fontsize='small')
            pl.savefig(plotfile)

        # Now accumulate pixel coordinates of corners of all new streaks to mask
        added_streak_ccds = []
        added_streak_corners = []

        for id, ccd in enumerate(new_ccdnum):
            ccd = new_ccdnum[id]
            if not maskit[id]:
                continue  # Only proceed with the ones to be masked
            # Get a pixel scale from the WCS, in arcsec/pix
            xmid = np.mean(ccd_corners_xpix)
            ymid = np.mean(ccd_corners_ypix)
            ra, dec = wcs[ccd].image2sky(xmid, ymid)
            ra2, dec2 = wcs[ccd].image2sky(xmid + 1, ymid)
            pixscale = np.hypot(np.cos(dec * np.pi / 180.) * (ra - ra2), dec - dec2)

            # width of streak, in pixels
            w = new_width[id] / pixscale + add_width
            if w <= 0.:
                continue  # Don't mask streaks of zero width
            # Make RA/Dec of track endpoints
            x = np.array([new_xc1[id], new_xc2[id]])
            y = np.array([new_yc1[id], new_yc2[id]])
            ra, dec = gnomonicInverse(x, y, ra0, dec0)
            # Convert to pixel coordinates
            x, y = wcs[ccd].sky2image(ra, dec)
            line = Line(x[0], y[0], x[1], y[1])
            # Create bounding rectangle of track
            corners_pix = boxTrack(line, w, ccd_xmin, ccd_xmax, ccd_ymin, ccd_ymax)
            added_streak_ccds.append(ccd)
            added_streak_corners.append(np.array(corners_pix))

        added_streak_ccds = np.array(added_streak_ccds)

        # Make new copies of streak files, adding new ones
        logger.debug('Rewriting streak files')

        for ccd, streakfile_in in streak_names.items():
            nmatch = len(re.findall(streak_name_in, streakfile_in))
            if nmatch != 1:
                logger.error('Could not update streak file named <' + streakfile_in + '>')
                return 1
            streakfile_out = re.sub(streak_name_in, streak_name_out, streakfile_in)
            # Use file system to make fresh copy of table's FITS file
            shutil.copy2(streakfile_in, streakfile_out)

            # Find new streaks for this ccd
            add_ids = np.where(added_streak_ccds == ccd)[0]
            if add_ids:
                # Open the table and add new streaks' info
                try:
                    fits = fitsio.FITS(streakfile_out, 'rw')
                    addit = np.recarray(len(add_ids),
                                        dtype=[('LABEL', '>i4'),
                                               ('CORNERS', '>f8', (4, 2)),
                                               ('CORNERS_WCS', '>f8', (4, 2))])
                    if fits[1]['LABEL'][:]:
                        first_label = np.max(fits[1]['LABEL'][:]) + 1
                    else:
                        first_label = 1
                    addit.LABEL = np.arange(first_label, first_label + len(addit))

                    for i, id in enumerate(add_ids):
                        corners_pix = added_streak_corners[id]
                        addit.CORNERS[i] = corners_pix
                        ra, dec = wcs[ccd].image2sky(corners_pix[:, 0], corners_pix[:, 1])
                        addit.CORNERS_WCS[i] = np.vstack((ra, dec)).T

                    fits[1].append(addit)
                    fits.close()
                except Exception as e:
                    print(e)
                    logger.error('Failure updating streak file <{:s}>'.format(streakfile_out))
                    return 1

        logger.debug('Remasking images')

        for imgfile_in in image_list:
            # Make the name needed for output
            nmatch = len(re.findall(image_name_in, imgfile_in))
            if nmatch != 1:
                logger.error('Could not create output name for image file named <' + imgfile_in + '>')
                return 1
            imgfile_out = re.sub(image_name_in, image_name_out, imgfile_in)

            logger.info(f"Loading image: {imgfile_in}")
            sci = DESImage.load(imgfile_in)
            ccd = sci.header['CCDNUM']

            # Find added streaks for this ccd
            add_ids = np.where(added_streak_ccds == ccd)[0]
            if add_ids:
                shape = sci.mask.shape
                yy, xx = np.indices(shape)
                points = np.vstack((xx.flatten(), yy.flatten())).T
                inside = None

                for id in add_ids:
                    # From Alex's immask routine: mark interior pixels
                    # for each added streak
                    v = added_streak_corners[id]
                    vertices = [(v[0, 0], v[0, 1]),
                                (v[1, 0], v[1, 1]),
                                (v[2, 0], v[2, 1]),
                                (v[3, 0], v[3, 1]),
                                (v[0, 0], v[0, 1])]
                    path = matplotlib.path.Path(vertices)

                    if inside is None:
                        inside = path.contains_points(points)
                    else:
                        inside = np.logical_or(inside, path.contains_points(points))

                # Make the list of masked pixels
                if inside is None:
                    ymask, xmask = np.array(0, dtype=int), np.array(0, dtype=int)
                else:
                    ymask, xmask = np.nonzero(inside.reshape(shape))

                sci.mask[ymask, xmask] |= parse_badpix_mask('STREAK')

            # Write something into the image header

            sci['DESCNCTS'] = time.asctime(time.localtime()) + \
                            ' Mask {:d} new streaks'.format(len(add_ids))
#            sci['HISTORY'] = time.asctime(time.localtime()) + \
#                             ' Mask {:d} new streaks'.format(len(add_ids))
            logger.info(f"Saving to: {imgfile_out}")
            sci.save(imgfile_out)

        logger.info('Finished connecting streaks')
        ret_code = 0
        return ret_code

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to streak connecting
        """
        parser.add_argument('--streak_file', type=str,
                            help='File holding list of input streak file names')
        parser.add_argument('--image_file', type=str,
                            help='File holding list of input image file names')
        parser.add_argument('--streak_name_in', type=str, default=cls.DEFAULT_STREAK_NAME_IN,
                            help='String to replace in input streak filenames')
        parser.add_argument('--streak_name_out', type=str, default=cls.DEFAULT_STREAK_NAME_OUT,
                            help='Replacement string for output streak filenames')
        parser.add_argument('--image_name_in', type=str, default=cls.DEFAULT_IMAGE_NAME_IN,
                            help='String to replace in input image filenames')
        parser.add_argument('--image_name_out', type=str, default=cls.DEFAULT_IMAGE_NAME_OUT,
                            help='Replacement string for output image filenames')
        parser.add_argument('--add_width', type=float, default=cls.DEFAULT_ADD_WIDTH,
                            help='number of pixels to grow (or shrink) streak width')
        parser.add_argument('--max_extrapolate', type=float, default=cls.DEFAULT_MAX_EXTRAPOLATE,
                            help='farthest to start a new streak from endpoint of an existing one (degrees)')
        parser.add_argument('--plotfile', type=str,
                            help='filename for diagnostic plot, if desired')

    @classmethod
    def run(cls, config):
        """Customized execution for streak connection.  No single input or output images.

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        streak_name_in = config.get(cls.step_name, 'streak_name_in')
        streak_name_out = config.get(cls.step_name, 'streak_name_out')
        image_name_in = config.get(cls.step_name, 'image_name_in')
        image_name_out = config.get(cls.step_name, 'image_name_out')
        add_width = config.getfloat(cls.step_name, 'add_width')
        max_extrapolate = config.getfloat(cls.step_name, 'max_extrapolate')

        if config.has_option(cls.step_name, 'plotfile'):
            plotfile = config.get(cls.step_name, 'plotfile')
        else:
            plotfile = None

        try:
            streak_list = filelist_to_list(config.get(cls.step_name, 'streak_file'))
        except:
            logger.error('Failure reading streak file names from {:s}'.format(streak_list))
            return 1

        try:
            image_list = filelist_to_list(config.get(cls.step_name, 'image_file'))
        except:
            logger.error('Failure reading image file names from {:s}'.format(image_list))
            return 1

        ret_code = cls.__call__(streak_list=streak_list,
                                image_list=image_list,
                                streak_name_in=streak_name_in,
                                streak_name_out=streak_name_out,
                                image_name_in=image_name_in,
                                image_name_out=image_name_out,
                                add_width=add_width,
                                max_extrapolate=max_extrapolate,
                                plotfile=plotfile)
        return ret_code


connect_streaks = ConnectStreaks()

# internal functions & classes

if __name__ == '__main__':
    connect_streaks.main()
