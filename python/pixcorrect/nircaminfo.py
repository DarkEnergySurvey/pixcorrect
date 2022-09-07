#!/usr/bin/env python3
"""
Collection of minutia supporting NIR camera data conversion to DECam-like.
Currently one instrument ESO VISTA...
"""

nir_paw_primary_keep=[
    'DATE',     
    'TELESCOP',
    'INSTRUME',
    'OBJECT',
    'RA',
    'DEC',
    'UTC',
    'LST',
    'PI-COI',
    'OBSERVER',
    'JITTR_ID',
    'NJITTER',
    'NOFFSETS',
    'OFFST_ID',
    'RECIPE',
    'REQTIME',
    'VSA_TIME',
    'VSA_MFID',
    'CASUVERS',
    'ESOGRADE',
    'OBSTATUS',
    'FILTER',
    'BAND',
    'CAMSYM',
    'EXPNUM']

nir_paw_hdu_keep=[
    'ORIGIN',
    'MJD-OBS',
    'DATE-OBS',
    'EXPTIME',
    'CTYPE1',
    'CTYPE2',
    'CRVAL1',
    'CRVAL2',
    'CRPIX1',
    'CRPIX2',
    'ORIGFILE',
    'CD1_1',
    'CD1_2',
    'CD2_1',
    'CD2_2',
    'PV2_1',
    'PV2_2',
    'PV2_3',
    'PV2_4',
    'PV2_5',
    'DARKCOR',
    'FLATCOR',
    'GAINCOR',
    'LINCOR',
    'CIR_CPM',
    'NDITCOR',
    'SKYLEVEL',
    'SKYNOISE',
    'CIR_VERS',
    'ELLIPTIC',
    'NUMBRMS',
    'STDCRMS',
    'WCSPASS',
    'SKYSUB',
    'DESTRIPE',
    'STRPRMS',
    'DRIBBLE',
    'PROV0000',
    'PROV0001',
    'PROV0002',
    'PROV0003',
    'PROV0004',
    'PROV0005',
    'PROV0006', 
    'PROV0007', 
    'NICOMB',
    'SEEING',
    'NUMSTDS',
    'PERCORR',
    'MAGZPT',
    'MAGZRR',
    'EXTINCT',
    'NUMZPT',
    'NIGHTZPT',
    'NIGHTZRR',
    'NIGHTNUM'
]

nir_vega_to_ab={
    'VZ':0.502,
    'VY':0.600,
    'J':0.916,
    'H':1.366,
    'Ks':1.827,
    'NB118':0.853
}

# This is set based on http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/linearity-sequences
nircam_satval={
    1: 33000.,
    2: 32000.,
    3: 33000.,
    4: 32000.,
    5: 24000.,
    6: 26000.,
    7: 35000.,
    8: 33000.,
    9: 35000.,
   10: 35000.,
   11: 37000.,
   12: 34000.,
   13: 33000.,
   14: 35000.,
   15: 34000.,
   16: 34000.
}


# https://manualzz.com/doc/o/a5c60/very-large-telescope-paranal-science-operations-vircam-6.3-filters (table 4)
# Note: McCracken etal 2012, A&A, 544, A156 (https://www.aanda.org/articles/aa/abs/2012/08/aa19507-12/aa19507-12.html)
#       gives ultraVista saturation as Y=14 and JHKs=15
nircam_satcalc={
    'VZ': 11.3,
    'VY': 10.8,
    'J': 11.1,
    'H': 11.0,
    'Ks': 10.2}


