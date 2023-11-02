#!/usr/bin/env python3
# Id: region_utils.py  
# CreatedBy: rgruendl  
# LastChangedBy: rgruendl 

"""Region handling utilities
"""


import time
import os
import re
import sys

#import matplotlib.path
#mport fitsio
import numpy as np


def loadRegionFile(region_file,reg_dict=None,AttChk=None,verbose=0):
    '''Read and parse region file into dict structure.
    Currently supports point and polygon type regions
    
    Inputs:
        region_file:   a DS9 region file expressing regions in sky coordinates
        reg_dict:      If exists then dictionary to be appended to (default=None)
        AttChk:        Dict of Attribute dictionaries to Check before allowing individual regions/lines into the dict
        verbose:       level ov verbosity (currently not used)

    Output:
        reg_dict       Dictionary of entries (one per region).

    Notes about region files: 
        Example formats for regions:
            fk5;pixel(ra,dec) # expnum=EXPNUM ccdnum=CCDNUM 
            fk5;polgon(ra1,dec1,ra2,dec2,ra3,dec3,....) # expnum=EXPNUM ccdnum=CCDNUM
        Key points:
            "fk5;" is not necessary... coordintes are currrently evaluated assuming they are in the same frame a the image
            "point" and "polygon" are used to distinguish how to further parse information.
             "(" and ")" denote the portion of the string to be parsed for coordinates
            ra,dec are assumed to be in units of degrees
            comments with format of expnum=VALUE ccdnum=VALUE are checked when AttChk is given 
    '''

    if (reg_dict is None):
        reg_dict={'nentry':0}

    print("Reading/parsing region mask file: {:s}".format(region_file))
    try:
        rfile=open(region_file,'r')
    except:
        print("Failed to read region file {:s}".format(region_file))
        exit(1)
    # Parse the region file        
    for line in rfile:
        if (line[0] != "#"):
            lbits=line.split()
#
#           If an attribute check is specified make sure that 
#           all attibutes are present and fulfill the check
#
            attchk_ok=True
            if (AttChk is not None):
                for att in AttChk:
                    fndmatch=False
                    attkey="{:s}=".format(att)
                    for bit in lbits:
                        if (re.match(attkey,bit)):
                            if (AttChk[att]['type']=='int'):
                                val=int(bit.split("=")[-1])
                            else:
                                val=bit.split("=")[-1]
                            if(AttChk[att]['val']==val):
                                fndmatch=True
                    if (not(fndmatch)):
                        attchk_ok=False
#
#           Attribute Check has passed.  Parse the line and add to reg_dict
#
            if (attchk_ok):
                n=reg_dict['nentry']+1
                reg_dict[n]={}
                reg_dict['nentry']=n
#                print("{:s}".format(line))
                c1=re.search("\(",line).end()
                c2=re.search("\)",line).start()
                if (re.search("point",line)):
                    reg_dict[n]['type']="point"
                if (re.search("polygon",line)):
                    reg_dict[n]['type']="polygon"
                radec=np.fromstring(line[c1:c2],dtype=float,sep=',')
                npts=int((radec.size)/2)
                reg_dict[n]['ra']=np.reshape(radec,(npts,2))[:,0]
                reg_dict[n]['dec']=np.reshape(radec,(npts,2))[:,1]
                reg_dict[n]['line']=line
#
#   Finished close file and return reg_dict
#
    rfile.close()          

    return reg_dict




