#! /bin/env python3

import argparse
import os
import re
import numpy as np
import fitsio
#import despyastro.wcsutil as wcsutil
import esutil.wcsutil as wcsutil
from pixcorrect import starmask_util as smu
import sys, copy
from collections import namedtuple
from scipy.spatial import ConvexHull



def get_bright_cat_objects(fname="coadd_objects.cat",wcs=None,config={'magcut':98.,'psf_rad':25.0}):
    cat_fits = fitsio.FITS(fname, 'r')
    objhdu=None
    for hname in ['OBJECTS','LDAC_OBJECTS']:
        if (hname in cat_fits):
            print("Found {:s}".format(hname))
            objhdu=hname
    if (objhdu is None):
        objhdu=1
    cat_cols=cat_fits[objhdu].get_colnames()
#    tdata=cat_fits[objhdu].read()
    cat_data = cat_fits[objhdu]['X_IMAGE','Y_IMAGE','MAG_AUTO','FLUX_RADIUS'][:]
    cat_fits.close()

    bright=np.where(np.logical_and(cat_data['MAG_AUTO']>0.0,cat_data['MAG_AUTO']<config['magcut']))
    x=cat_data['X_IMAGE'][bright]
    y=cat_data['Y_IMAGE'][bright]
    r=np.minimum(4.0*cat_data['FLUX_RADIUS'][bright],config['psf_rad'])
#    r=4.0*cat_data['FLUX_RADIUS'][bright]
    print("After magnitude cut there are {:d} OBJECTS to mask".format(x.size))
    if (wcs is None):
        for i in range(x.size):
            print("image;circle({:.2f},{:.2f},{:.2f}) # color=cyan width=3 BrightObj ".format(x[i],y[i],r[i]))
    else: 
        ra,dec=wcs.image2sky(x,y)
        for i in range(ra.size):
            print("fk5;circle({:.7f},{:.7f},{:.2f}\") # color=cyan width=3 BrightObj ".format(ra[i],dec[i],r[i]*0.263))

    return [x,y,r]


def excess_noise_cut(g_mag,slope=0.3,intercept=0.15):
     """ Cut on astrometric_excess_noise. Returns log10(astrometric_excess_noise) """ 
     log_noise = slope*(g_mag - 18.2) + intercept
     return np.clip(log_noise,0.3,None)


def get_bright_DB_objects(dbTable,dbSection='db-dessci',dbSchema='des_admin',coadd_head=None,wcs=None,config={'psf_rad':25.0},verbose=0):

    import despydb.desdbi

    radec_box={}
    if ('CROSSRA0' in coadd_head):
        if (coadd_head['CROSSRA0']=="Y"):
            radec_box['crossra0']=True
        else:
            radec_box['crossra0']=False
        radec_box['ra1']=coadd_head['RACMIN']
        radec_box['ra2']=coadd_head['RACMAX']
        radec_box['dec1']=coadd_head['DECCMIN']
        radec_box['dec2']=coadd_head['DECCMAX']
    else:
        if ('NAXIS1' in coadd_head):
            nx=coadd_head['NAXIS1']
            ny=coadd_head['NAXIS2']
        else:
            print("No COADD header information?  Aborting!")
            exit(1)
        x=np.array([1,1,nx,nx])
        y=np.array([1,ny,ny,1])
        if (wcs is None):
            print("No COADD WCS information?  Aborting!")
            exit(1)
        ra,dec=wcs.image2sky(x,y)
        radec_box['ra1']=ra.amax
        radec_box['ra2']=ra.amin
        if (radec_box['ra1']-radec_box['ra2']>180.0):
            radec_box['crossra0']=True
        else:
            radec_box['crossra0']=False
        radec_box['dec1']=dec.amin
        radec_box['dec2']=dec.amax

    radec_box=smu.expand_range(radec_box,method='fixed',extend=0.5,verbose=verbose)

#
#   Prepare a database connection to obtain catalog info
#   Make an RA-Dec box query of 2MASS and get positions with J-band magnitudes.
#
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile, dbSection, retry=True)

#    if (use2MASS):
#        print("Querying 2MASS_PSC")
#        StarCat, StarHead=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='TWOMASS_PSC',cols=['RA','DEC','J'],Timing=False,verbose=0)
#        CatMag='J'
##    else:
    if (dbTable == "GAIA_DR2"):
        print("Querying GAIA_DR2")
        StarCat, StarHead=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='GAIA_DR2',cols=['RA','DEC','PHOT_G_MEAN_MAG'],Timing=False,verbose=0)
        CatMag='PHOT_G_MEAN_MAG'
    elif (dbTable == "GAIA_EDR3"):
        print("Querying GAIA_EDR3")
        StarCat, StarHead=smu.get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='GAIA_EDR3',cols=['RA','DEC','PHOT_G_MEAN_MAG','ASTROMETRIC_EXCESS_NOISE'],Timing=False,verbose=0)
        CatMag='PHOT_G_MEAN_MAG'
#
#       Remove galaxies
#
        print(StarHead)
        gaia_stars = np.where(np.log10(StarCat['ASTROMETRIC_EXCESS_NOISE']) < excess_noise_cut(StarCat['PHOT_G_MEAN_MAG']))
        print(" Pre-galaxy cut: {:d}".format(StarCat['RA'].size))
        for col in StarHead:
            StarCat[col]=StarCat[col][gaia_stars]
        print("Post-galaxy cut: {:d}".format(StarCat['RA'].size))



    if (wcs is not None):
        x,y=wcs.sky2image(StarCat['RA'],StarCat['DEC'])
        r=config['psf_rad']*np.ones(x.shape)
    else:
        r=config['psf_rad']*np.ones(StarCat['RA'].shape)

    if (wcs is None):
        for i in range(StarCat['RA'].size):
            print("fk5;circle({:.7f},{:.7f},{:.2f}\") # color=cadet blue width=3 BrightObj ".format(StarCat['RA'][i],StarCat['DEC'][i],r[i]*0.263))
    else: 
        for i in range(x.size):
            print("fk5;circle({:.7f},{:.7f},{:.2f}\") # color=cadet blue width=3 BrightObj ".format(StarCat['RA'][i],StarCat['DEC'][i],r[i]*0.263))

    return [x,y,r]



def delete_pix(x,y,r,pix):
#   Remove pixels that "collide" with given circular region (e.g. a region around a bright star)
    for i in range(len(x)):
        del_px = np.where((x[i]-pix[:,1])**2+(y[i]-pix[:,2])**2<r[i]**2)[0]
        pix=np.delete(pix,del_px,0)
    return pix


# https://stackoverflow.com/questions/62980280/finding-neighboring-pixels-python
def candidate_neighbors(node):
    return ((node[0]-1, node[1]-1), (node[0]-1, node[1]), (node[0]-1, node[1]+1), (node[0], node[1]-1), 
            (node[0], node[1]+1), (node[0]+1, node[1]-1), (node[0]+1, node[1]), (node[0]+1, node[1]+1))


def neighboring_groups(nodes):
    remain = set(nodes)
    while len(remain) > 0:
        visit = [remain.pop()]
        group = []
        while len(visit) > 0:
            node = visit.pop()
            group.append(node)
            for nb in candidate_neighbors(node):
                if nb in remain:
                    remain.remove(nb)
                    visit.append(nb)
        yield tuple(group)


def form_output_head(head):
    # form a minimal output header (tested for TAN and maybe ready for TPV projections)

    std_key=['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CD1_1','CD1_2','CD2_1','CD2_2']
    tpv_key=['PV1_0','PV1_1','PV1_2','PV1_3','PV1_4','PV1_5','PV1_6','PV1_7','PV1_8','PV1_9','PV1_10',
             'PV2_0','PV2_1','PV2_2','PV2_3','PV2_4','PV2_5','PV2_6','PV2_7','PV2_8','PV2_9','PV2_10']
    out_head={}
    isTPV=False
    if ('CTYPE1' in head):
        if (re.match("TPV",head['CTYPE1'])):
            isTPV=True
    for hitem in head:
        if (hitem in std_key):
            out_head[hitem]=head[hitem]
        if ((isTPV)and(hitem in tpv_key)):
            out_head[hitem]=head[hitem]

    return out_head

    
def expand_hull(hull,hull_cen=None,exp=0.):
    # Expand a hull by factor 1+exp.
    # If center is not give then use average from vertices (better to get it from points that generated hull)
    if (hull_cen is None):
        hull_cen=np.mean(hull,axis=0)
    exp_hull=hull_cen+((1.+exp)*(hull-hull_cen))
    return exp_hull    
        

def expand_pix(pix, exp, dims, wcs=None, frame=-1):
    # Currently deprecated (in favor of group_to_region and then polygon_to_pix below...)
    # grow the mask region
    if (frame is None):
        fnum=-1
    else:
        fnum=frame
    out_pix = np.flip(np.array(pix, dtype=int), 1)
    print("Begin work on pixel group with {:d} pixels".format(out_pix.shape[0]))

    # RegionDebug
    if (config['debug']):
        print("RAG expand: ",out_pix[:,1],out_pix[:,0])
        ra,dec=wcs.image2sky(out_pix[:,1],out_pix[:,0])
        for i in range(ra.size):
            print("fk5;point({:.7f},{:.7f}) # point=x color=green EFrame {:d} Z".format(ra[i],dec[i],fnum))
#
#   Try to fit a hull (so that a polygon can be defined for the region.
#   Uses:
#       https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
#       to get the convex hull
#       https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
#

    try:
        hull = out_pix[ConvexHull(out_pix).vertices,:]
        if ((config['debug'])or(config['verbose']>2)):
            print("Found hull: ",hull)
        else:
            if (config['verbose']>0):
                print("Found hull")
    except Exception as e:
        print("Qhull Error Encountered: {:}".format(e))
        print("pixel list: ",out_pix)       
        ra,dec=wcs.image2sky(out_pix[:,1],out_pix[:,0])
        # RegionKeep
        for i in range(ra.size):
            print("fk5;point({:.7f},{:.7f}) # point=x color=red EFrame {:d} Z".format(ra[i],dec[i],fnum))
        return out_pix

#   Hull succeeded:

    pts_sz=np.max(out_pix,axis=0)-np.min(out_pix,axis=0)
    pts_cen=np.mean(out_pix,axis=0)
    if ((config['debug'])or(config['verbose']>2)):
        print("pts_cen: ",pts_cen)
    # RAG RegionDebug
    if (config['debug']):
        ra,dec=wcs.image2sky(hull[:,1],hull[:,0])
        radec_vert=["{:.7f},{:.7f}".format(ra[i],dec[i]) for i in range(ra.size)]
        print("fk5;polygon({:s}) # color=green EFrame {:d} Z".format(",".join(radec_vert),fnum))

    exp_hull=expand_hull(hull,hull_cen=pts_cen,exp=exp)
    # RAG RegionKeep
    ra,dec=wcs.image2sky(exp_hull[:,1],exp_hull[:,0])
    radec_vert=["{:.7f},{:.7f}".format(ra[i],dec[i]) for i in range(ra.size)]
    print("fk5;polygon({:s}) # color=orange width=2 EFrame {:d} Z".format(",".join(radec_vert),fnum))
        
    # define the set of enclosed pixels. - do I need a 1px offset here????
    ####### CHECK ############
    exp_hull_min=np.rint(np.amin(exp_hull,axis=0))-2
    exp_hull_max=np.rint(np.amax(exp_hull,axis=0))+2
    if (exp_hull_min[0]<0):
        exp_hull_min[0]=0
    if (exp_hull_min[1]<0):
        exp_hull_min[1]=0
    if (exp_hull_max[0]>dims[0]):
        exp_hull_max[0]=dims[0]
    if (exp_hull_max[1]>dims[1]):
        exp_hull_max[1]=dims[1]

#   This is a quick guess (may need to be fixed)

    xy=polygon_to_pix(exp_hull[:,1],exp_hull[:,0])

#   Also the outputs can now exceed an image boundary so need to fix that too

    return xy
 

def polygon_to_pix(x,y):
    # From a set of x,y coordinates expressing vertices of polygon
    # Determine the set of x,y pixels within that region.
    #
    # NOTE: The current implementation assumes the vertices describe
    #       a convex polygon (hull).

    xvt1=np.zeros(x.size)
    xvt2=np.zeros(x.size)
    yvt1=np.zeros(x.size)
    yvt2=np.zeros(x.size)
    mside=np.zeros(x.size)
    bside=np.zeros(x.size)
    nvt=x.size

#   First define set of lines (vertices, slopes, and intercepts) whose segments bound the hull.

    for i in range(nvt):
        xvt1[i]=x[i]
        yvt1[i]=y[i]
        if (i == nvt-1):
            xvt2[i]=x[0]
            yvt2[i]=y[0]
        else:
            xvt2[i]=x[i+1]
            yvt2[i]=y[i+1]

        if (xvt2[i]!=xvt1[i]):
            mside[i]=(yvt2[i]-yvt1[i])/(xvt2[i]-xvt1[i])
            bside[i]=((yvt1[i]*xvt2[i])-(yvt2[i]*xvt1[i]))/(xvt2[i]-xvt1[i])
        else:
            mside[i]=None
            bside[i]=None
#
#   Now take the vertex endpoints and make them strictly increasing (so they can now be re-used as ranges)
#
    for i in range(nvt):
        if (xvt1[i]>xvt2[i]):
            xtmp=xvt1[i]
            xvt1[i]=xvt2[i]
            xvt2[i]=xtmp
        if (yvt1[i]>yvt2[i]):
            ytmp=yvt1[i]
            yvt1[i]=yvt2[i]
            yvt2[i]=ytmp
#        if (mside[i] == 0.0):
#            print("Wall",i," at ",yvt1[i],yvt2[i])

    # define the set of enclosed pixels. - do I need a 1px offset here????
    ####### CHECK ############
#    exp_hull_min=np.rint(np.amin(exp_hull,axis=0))-2
#    exp_hull_max=np.rint(np.amax(exp_hull,axis=0))+2
#    if (exp_hull_min[0]<0):
#        exp_hull_min[0]=0
#    if (exp_hull_min[1]<0):
#        exp_hull_min[1]=0
#    if (exp_hull_max[0]>dims[0]):
#        exp_hull_max[0]=dims[0]
#    if (exp_hull_max[1]>dims[1]):
#        exp_hull_max[1]=dims[1]
    xmin=np.rint(np.amin(x))-2
    xmax=np.rint(np.amax(x))+2
#    ymin=np.rint(np.amin(y))-2
#    ymax=np.rint(np.amax(y))+2

#   Finally step through the range covered by polygon/hull and flag the points      
    mask_px=[]
    for ix in np.arange(xmin,xmax+1):
#        print("Starting ix=",ix)
        yint=[]
        for i in range(nvt):
            if (np.isnan(mside[i])):
                if (ix == np.rint(xvt1[i])):
#                    print("Y-range:",yvt1[i],yvt2[i])
                    yint.append(yvt1[i])
                    yint.append(yvt2[i])
            elif (mside[i] == 0):
                if ((ix > xvt1[i])and(ix <= xvt2[i])):
                    yint.append(yvt2[i])
            else:
                yval=mside[i]*ix+bside[i]
                if ((yval >= yvt1[i])and(yval <= yvt2[i])):
#                    print("Hit ",i," at ",yval)
                    yint.append(yval)
        yint=np.array(yint)
        if (yint.size > 0):
#            print(" {:.2f}  {:.2f} {:.2f}  {:d} ".format(ix,np.amin(yint),np.amax(yint),yint.size))
#            yset=[]
            ymin=np.rint(np.amin(yint))
            ymax=np.rint(np.amax(yint))+1
#            if (ymin < 0):
#                ymin=0
##            if (ymax > dims[1]-1):
##                ymax=dims[1]-1
##            for iy in np.arange(np.rint(np.amin(yint)),np.rint(np.amax(yint))+1):
#            if (ymax > dims[1]):
#                ymax=dims[1]
            for iy in np.arange(ymin,ymax):
#                yset.append(iy)
                mask_px.append((ix,iy))
#            print(yset)
#        else:
#            print(" {:.2f}  No points ".format(ix))

#
#   There is a rare corner case for a small polygon where no pixels will be found. 
#   In this case flag the center pixel and the eight surrounding pixels.
#
    if (len(mask_px)<1):
        print("Small polygon corner case detected.... treating as a point using the centroid of the vertices")
        xpt_exp=np.array([1,0,-1,1,0,-1,1,0,-1])
        ypt_exp=np.array([1,1,1,0,0,0,-1,-1,-1])
        xmsk=np.rint(np.mean(xvt1))+xpt_exp
        ymsk=np.rint(np.mean(yvt1))+ypt_exp
        mask_px=[(xmsk[i],ymsk[i]) for i in range(xmsk.size)]
    
    xy = np.array([(pt[0], pt[1]) for pt in mask_px])

    return np.array(xy, dtype=int)


def group_to_region(pix, exp, dims, wcs=None ,frame=-1, debug=False, verbose=0):
    # Take a pixel group
    # Find the convex hull that encompasses it (if it fails to find a hull then write individual points)
    # Expand hull if requested
    # Place results in dict for external consumption.

    if (frame is None):
        fnum=-1
    else:
        fnum=frame

    rdict={}
    rdict['nentry']=0
    rdict['npoint']=0
    rdict['npoly']=0
    rdict['area']=0

    out_pix = np.flip(np.array(pix, dtype=int), 1)
    print("Begin work on pixel group with {:d} pixels".format(out_pix.shape[0]))

    # RegionDebug
    if (debug):
        print("Workiing with group of pixels: ",out_pix[:,1],out_pix[:,0])
        ra,dec=wcs.image2sky(out_pix[:,1],out_pix[:,0])
        for i in range(ra.size):
            print("fk5;point({:.7f},{:.7f}) # point=x color=green EFrame {:d} Z".format(ra[i],dec[i],fnum))
#
#   Try to fit a hull (so that a polygon can be defined for the region.
#   Uses:
#       https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
#       to get the convex hull
#       https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
#
    try:
        hull = out_pix[ConvexHull(out_pix).vertices,:]
        if ((debug)or(verbose>2)):
            print("Found hull: ",hull)
    except Exception as e:
        print("Qhull Error Encountered: {:}".format(e))
        print("pixel list: ",out_pix)       
        ra,dec=wcs.image2sky(out_pix[:,1],out_pix[:,0])
        # RegionKeep
        for i in range(ra.size):
            print("fk5;point({:.7f},{:.7f}) # point=x color=red EFrame {:d} Z".format(ra[i],dec[i],fnum))
            rdict['nentry']=rdict['nentry']+1
            rdict['npoint']=rdict['npoint']+1
            rdict[rdict['nentry']]='fk5;point({:.7f},{:.7f}) # point=x'.format(ra[i],dec[i])
        return rdict

#   Hull succeeded:

    pts_sz=np.max(out_pix,axis=0)-np.min(out_pix,axis=0)
    pts_cen=np.mean(out_pix,axis=0)
    if ((debug)or(verbose>2)):
        print("pts_cen: ",pts_cen)
    # RAG RegionDebug
    if (debug):
        ra,dec=wcs.image2sky(hull[:,1],hull[:,0])
        radec_vert=["{:.7f},{:.7f}".format(ra[i],dec[i]) for i in range(ra.size)]
        print("fk5;polygon({:s}) # color=green EFrame {:d} Z".format(",".join(radec_vert),fnum))

    if (exp == 0):
        exp_hull=hull
    else:
        exp_hull=expand_hull(hull,hull_cen=pts_cen,exp=exp)
    # RAG RegionKeep
    ra,dec=wcs.image2sky(exp_hull[:,1],exp_hull[:,0])
    radec_vert=["{:.7f},{:.7f}".format(ra[i],dec[i]) for i in range(ra.size)]
    print("fk5;polygon({:s}) # color=orange width=2 EFrame {:d} Z".format(",".join(radec_vert),fnum))
    rdict['nentry']=rdict['nentry']+1
    rdict['npoly']=rdict['npoly']+1
    rdict['area']=rdict['area']+out_pix.shape[0]
    rdict[rdict['nentry']]='fk5;polygon({:s}) # '.format(",".join(radec_vert))

    return rdict

