# $Id: starmask_utils.py 46990 2018-05-10 19:57:25Z rgruendl $
# $Rev:: 46990                            $:  # Revision of last commit.
# $LastChangedBy:: rgruendl               $:  # Author of last commit.
# $LastChangedDate:: 2018-05-10 14:57:25 #$:  # Date of last commit.

"""Starmask utility functions
"""

import re
import pandas as pd
import time
import numpy as np
from pixcorrect.corr_util import logger
import despydb.desdbi

###########################################
def expand_range(radec_box,method='fixed',extend=10.,verbose=0):

    """Expand the size of an radec_box, either by a fixed or fractional amount
       Check to make sure that the extent has not crossed the RA=0h divide (if so update appropriately)
    """
#
#       Stolen code from mepipelineappintg.cat_query (left it close to the same so that it could be integrated together by
#       moving functionality somewhere else.
#
    if (method == 'fractional'):
        print("Expansion method: fractional ({:.2f} percent)".format(100.*extend))
        if (radec_box['crossra0']):
            dra=radec_box['ra2']-(radec_box['ra1']-360.0)
        else:
            dra=radec_box['ra2']-radec_box['ra1']
        ddec=radec_box['dec2']-radec_box['dec1']

        rmin=radec_box['ra1']-(extend*dra)
        rmax=radec_box['ra2']+(extend*dra)
        dmin=radec_box['dec1']-(extend*ddec)
        dmax=radec_box['dec2']+(extend*ddec)

    elif (method == 'fixed'):
        print("Expansion method: fixed ({:.1f} arcminutes)".format(extend))
        deccen=0.5*(radec_box['dec2']+radec_box['dec1'])
        cosdec=np.cos(deccen*np.pi/180.)
        dra=extend/60.0/cosdec
        ddec=extend/60.0

        rmin=radec_box['ra1']-dra
        rmax=radec_box['ra2']+dra
        dmin=radec_box['dec1']-ddec
        dmax=radec_box['dec2']+ddec

    else:
        print("Warning! Unrecognized method: '{:s}'.  No changes made to search range.".format(method))
        rmin=radec_box['ra1']
        rmax=radec_box['ra2']
        dmin=radec_box['dec1']
        dmax=radec_box['dec2']

#
#   Now some simple checking
#
    if (verbose > 0):
        print("Expanded: {:}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
            radec_box['crossra0'],
            rmin,radec_box['ra1'],radec_box['ra2'],rmax,
            dmin,radec_box['dec1'],radec_box['dec2'],dmax))

#
#   In the unlikely even that the North/South Celestial Pole was crossed just truncate it at the pole (do not deal with the unholy case that really occurred)
#
    if (dmin<-90.):
        dmin=-90.0
        if (verbose > 0):
            print("Ammended: {:}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                radec_box['crossra0'],
                rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                dmin,radec_box['dec1'],radec_box['dec2'],dmax))
    if (dmax>90.):
        dmax=90.0
        if (verbose > 0):
            print("Ammended: {:}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                radec_box['crossra0'],
                rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                dmin,radec_box['dec1'],radec_box['dec2'],dmax))

#
#   Check to make sure that if RA range crosses the RA=0/24 boundary is handled consistent with COADDTILE_GEOM/IMAGE table structure
#
    if (not(radec_box['crossra0'])):
        if (rmin < 0.):
            rmin=rmin+360.0
            radec_box['crossra0']=True
            if (verbose > 0):
                print("Ammended: {:}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                    radec_box['crossra0'],
                    rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                    dmin,radec_box['dec1'],radec_box['dec2'],dmax))
        if (rmax > 360.):
            rmax=rmax-360.0
            radec_box['crossra0']=True
            if (verbose > 0):
                print("Ammended: {:}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                    radec_box['crossra0'],
                    rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                    dmin,radec_box['dec1'],radec_box['dec2'],dmax))

    radec_box['ra1']=rmin
    radec_box['ra2']=rmax
    radec_box['dec1']=dmin
    radec_box['dec2']=dmax
#
#   Check to make sure that if RA range crosses the RA=0/24 boundary is handled consistent with COADDTILE_GEOM/IMAGE table structure
#
    if (not(radec_box['crossra0'])):
        if (rmin < 0.):
            rmin=rmin+360.0
            radec_box['crossra0']=True
            if (verbose > 0):
                print("Ammended: {:s}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                    radec_box['crossra0'],
                    rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                    dmin,radec_box['dec1'],radec_box['dec2'],dmax))
        if (rmax > 360.):
            rmax=rmax-360.0
            radec_box['crossra0']=True
            if (verbose > 0):
                print("Ammended: {:s}  RA: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f}   Dec: {:9.5f} << {:9.5f} -- {:9.5f} >> {:9.5f} ".format(
                    radec_box['crossra0'],
                    rmin,radec_box['ra1'],radec_box['ra2'],rmax,
                    dmin,radec_box['dec1'],radec_box['dec2'],dmax))
#
#
#
    radec_box['ra1']=rmin
    radec_box['ra2']=rmax
    radec_box['dec1']=dmin
    radec_box['dec2']=dmax

    return radec_box

###########################################
def get_cat_radec_range(radec_box,dbh,dbSchema='des_admin.',table='GAIA_DR2',cols=['ra','dec','phot_g_mean_mag'],Timing=False,verbose=0):

    """ Pull Catalog Data in an RA/Dec range. Default is GAIA_DR2 objects (ra,dec,mag) in a region

    radec_box:  Dict w/ keys: ra1,ra2,dec1,dec2, and crossra0[bool] that describe box to search
    dbh:        Database connection to be used
    dbSchema:   DB schema (default='DES_ADMIN')
    table:      Catalog Table (must have RA,Dec) to query for objects (default='GAIA_DR2')
    cols:       List of columns to return (default is ['ra','dec','phot_g_mean_mag'])
    verbose:    Sets level of verbosity."

    Return:     Dict of numpy arrays (one for each column), List of Columns

    """

    t0=time.time()
    if (radec_box['crossra0']):
#
#       Form Query for case where RA ranges crosses RA=0h (not very good at poles)
#
        query="""select {cname:s}
            from {schema:s}{tbl:s}
            where (ra < {r2:.6f} or ra > {r1:.6f})
                and dec between {d1:.6f} and {d2:.6f}""".format(
        cname=",".join(cols),
        schema=dbSchema,
        tbl=table,
        r1=radec_box['ra1'],
        r2=radec_box['ra2'],
        d1=radec_box['dec1'],
        d2=radec_box['dec2'])
    else:
#
#       Form query for normal workhorse case
#
        query="""select {cname:s}
            from {schema:s}{tbl:s}
            where ra between {r1:.6f} and {r2:.6f}
                and dec between {d1:.6f} and {d2:.6f}""".format(
        cname=",".join(cols),
        schema=dbSchema,
        tbl=table,
        r1=radec_box['ra1'],
        r2=radec_box['ra2'],
        d1=radec_box['dec1'],
        d2=radec_box['dec2'])
#
    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print("{:s}".format(QueryOneLine))
        if (verbose > 1):
            print("{:s}".format(query))
#
#   Establish a DB cursor
#
    curDB = dbh.cursor()

    prefetch=100000
    curDB.arraysize=int(prefetch)
    curDB.execute(query)
    header=[d[0].upper() for d in curDB.description]
    cat_data=pd.DataFrame(curDB.fetchall())

    CatDict={}
    if (cat_data.empty):
        print("# No values returned from query of {tval:s} ".format(tval="GAIA_DR2"))
        for val in header:
            CatDict[val]=np.array([])
    else:
        cat_data.columns=header
        for val in header:
            CatDict[val]=np.array(cat_data[val])
    curDB.close()

    if (verbose>0):
        print("# Number of objects found in {schema:s}{tbl:s} is {nval:d} ".format(
            schema=dbSchema,
            tbl=table,
            nval=CatDict[header[0]].size))
    if (Timing):
        t1=time.time()
        print(" Query execution time: {:.2f}".format(t1-t0))

    return CatDict,header
