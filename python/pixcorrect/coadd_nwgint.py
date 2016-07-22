#!/usr/bin/env python

from pixcorrect.null_weights import null_weights
from pixcorrect.row_zipper   import row_zipper
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

from despyastro.CCD_corners import update_DESDM_corners

import despyfits
from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage,update_hdr_compression,insert_eupspipe
from despymisc.miscutils import elapsed_time
from despyfits import updateWCS 
import time

import fitsio
import numpy as np

class CoaddZipperInterpNullWeight(PixCorrectMultistep):

    """Run custom weights for STAR and do not write MSK plane for multi-epoch (me)'"""

    config_section = "coadd_nwgit"
    description = 'Perform zipper interpolation along rows and null_weights in one step'
    step_name = config_section
    DEFAULT_CUSTOM_WEIGHT = False
    DEFAULT_HEADFILE = False
    DEFAULT_HDUPCFG = False
    DEFAULT_TILENAME = False
    DEFAULT_TILEID = False
    DEFAULT_ME_WGT_KEEPMASK = False
    DEFAULT_NULL_MASK_SCI = '0'

    # Fix the step_name for passing the command-line arguments to the classes
    null_weights.__class__.step_name = config_section
    row_zipper.__class__.step_name   = config_section
    
    def __call__(self):
        """
        Run row_zipper and null_weights in one step, we run the tasks
        by calling step_run in each class
        """
        t0 = time.time()

        # Check if we want special multi-epoch weighting, and which bits we want to 'save'
        me_wgt_keepmask = get_safe_boolean('me_wgt_keepmask',self.config,self.config_section)

        # Get verbose
        try:
            verbose = self.config.get(self.config_section,'verbose')
        except:
            verbose = False

        # Get the science image
        input_image = self.config.get(self.config_section,'in')
        self.sci = DESImage.load(input_image)

        # Add TILENAME and TILEID to sci header (optional) if required
        self.update_sci_header(input_image)

        # Update the header wcs if both headfile and hdupcfg are present (optional)
        self.update_wcs_header(input_image,verbose=verbose)
        
        # Check if want to create the custon weight for SWArp/SExtractor combination
        if me_wgt_keepmask:
            self.custom_weight(input_image)
        
        # Run null_weights
        t1 = time.time()
        logger.info("Running null_weights on: %s" % input_image)
        null_weights.step_run(self.sci,self.config)
        logger.info("Time NullWeights : %s" % elapsed_time(t1))

        # Run row_zipper
        t2 = time.time()
        logger.info("Running row_zipper on: %s" % input_image)
        row_zipper.step_run(self.sci,self.config)
        logger.info("Time ZipperInterp : %s" % elapsed_time(t2))

        # Null the sci image only if null_mask_sci !=0
        self.null_sci(input_image)
        
        output_image = self.config.get(self.config_section, 'out')
        # Special write out
        if me_wgt_keepmask :
            self.custom_write(output_image)
        else:
            self.sci.save(output_image)
        
        logger.info("Wrote new file: %s" % output_image)
        logger.info("Time Total: %s" % elapsed_time(t0))

        return 0

    def update_wcs_header(cls,input_image,verbose=False):

        # Get optional config file, first we try to get them as boolean, then as strings
        headfile = get_safe_boolean('headfile',cls.config,cls.config_section)
        hdupcfg  = get_safe_boolean('hdupcfg',cls.config,cls.config_section)

        # Update the header if both headfile and hdupcfg are present
        if  headfile and hdupcfg:
            logger.info("Will update image header with scamp .head file %s" % headfile)
            cls.sci = updateWCS.run_update(cls.sci,headfile=headfile,hdupcfg=hdupcfg,verbose=verbose)
    
    def update_sci_header(cls,input_image):
        
        tilename = get_safe_boolean('tilename',cls.config,cls.config_section)
        tileid   = get_safe_boolean('tileid',cls.config,cls.config_section)
        if tilename:
            record={'name':'TILENAME', 'value':tilename, 'comment':'DES Tilename'}
            cls.sci.header.add_record(record)
        if tileid:
            record={'name':'TILEID', 'value':int(tileid), 'comment':'Tile ID for DES Tilename'}
            cls.sci.header.add_record(record)


    def null_sci(cls, input_image):
        
        null_mask_sci = parse_badpix_mask( cls.config.get(cls.config_section, 'null_mask_sci') )
        if null_mask_sci !=0:
            logger.info('Nulling science image from null_mask_bits')
            kill  = np.array(cls.sci.mask & null_mask_sci, dtype=bool)
            cls.sci.data[kill]  = 0.0
        else:
            logger.info('Science image was not null')

        return

    def custom_weight(cls,input_image):
        # Make custom weight, that will not zero STAR maskbit
        logger.info("Will perform special weighting for multi-epoch input on %s" % input_image)
        # Make a copy of the original untouched weight
        cls.sci.weight_custom = np.copy(cls.sci.weight)
        null_mask       = cls.config.get(cls.config_section, 'null_mask')
        me_wgt_keepmask = cls.config.get(cls.config_section, 'me_wgt_keepmask')

        # Make python lists of the coma-separated input lists
        null_list = null_mask.split(',')
        keep_list = me_wgt_keepmask.split(',') 

        # Special case we care:
        # . we are nulling the TRAIL but want keep where STAR 
        if 'TRAIL' in null_list and 'STAR' in keep_list and 'TRAIL' not in keep_list:
            # Remove STAR from the list
            if 'STAR' in null_list: null_list.remove('STAR')
            null_mask_bits = parse_badpix_mask(','.join(null_list))
            # Null each plane at a time. First the TRAILS and replace with STAR
            kill  = np.array(cls.sci.mask & parse_badpix_mask('TRAIL'), dtype=bool)
            stars = np.array(cls.sci.mask & parse_badpix_mask('STAR'), dtype=bool)
            cls.sci.weight_custom[kill]  = 0.0
            cls.sci.weight_custom[stars]  = np.copy(cls.sci.weight[stars])
            # Loop over the bitplanes, but skipping TRAIL, which we already did
            null_list.remove('TRAIL')
            for bitplane in null_list:
                kill  = np.array(cls.sci.mask & parse_badpix_mask(bitplane), dtype=bool)
                cls.sci.weight_custom[kill]  = 0.0
        # We  remove tham from the null_list
        else:
            for bitplane in me_wgt_keepmask.split(','):
                if bitplane in null_list: null_list.remove(bitplane)
            null_mask_bits = parse_badpix_mask(','.join(null_list))
            kill = np.array( cls.sci.mask & null_mask_bits, dtype=bool)
            cls.sci.weight_custom[kill]  = 0.0
            
    def custom_write(cls,output_image):
        # Write out the image using fitsio, but skipping the mask as we won't need it.
        ofits = fitsio.FITS(output_image,'rw',clobber=True)

        # Here we mimick the steps followed by DESImage.save()
        # SCI
        logger.info("Creating SCI HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        cls.sci.header = update_hdr_compression(cls.sci.header,'SCI')
        logger.info("Calculating CCD corners/center/extern keywords for SCI HDU ")
        cls.sci.header = update_DESDM_corners(cls.sci.header,get_extent=True, verb=False)
        if despyfits.DESImage.pipekeys_write:
            logger.info("Inserting EUPS PIPEPROD and PIPEVER to SCI HDU")
            cls.sci.header = insert_eupspipe(cls.sci.header)
        ofits.write(cls.sci.data,  extname='SCI', header=cls.sci.header)
        # WGT
        logger.info("Creating WGT HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        cls.sci.weight_hdr = update_hdr_compression(cls.sci.weight_hdr,'WGT')
        ofits.write(cls.sci.weight,extname='WGT',header=cls.sci.weight_hdr)
        # WGT_ME 
        # For  WGT_ME we do not need to update the FZ keywords, as we use the same hdr as WGT
        logger.info("Creating WGT_ME HDU")
        ofits.write(cls.sci.weight_custom,extname='WGT_ME',header=cls.sci.weight_hdr)
        # MSK
        logger.info("Creating MSK HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        cls.sci.mask_hdr = update_hdr_compression(cls.sci.mask_hdr,'MSK')
        ofits.write(cls.sci.mask,extname='MSK',header=cls.sci.mask_hdr)
        ofits.close()

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments for null_weights and row_zipper
        """
        null_weights.add_step_args(parser)
        row_zipper.add_step_args(parser)
        #parser.add_argument('--custom_weight', action='store_true',default=cls.DEFAULT_CUSTOM_WEIGHT,
        #                    help='Run custom weights for STAR and do not write MSK plane for multi-epoch (me)')
        parser.add_argument('--me_wgt_keepmask', action='store',default=cls.DEFAULT_ME_WGT_KEEPMASK,
                            help='Run custom weight for multi-epoch (me) WGT_ME and preserve KEEPMASK')
        parser.add_argument('--null_mask_sci', action='store',default=cls.DEFAULT_NULL_MASK_SCI,
                            help='Names of mask bits to null (or an integer mask) on the SCI plane')
        parser.add_argument('--headfile', action='store', default=cls.DEFAULT_HEADFILE,
                            help='Headfile (containing most update information)')
        parser.add_argument('--hdupcfg', action='store', default=cls.DEFAULT_HDUPCFG,
                            help='Configuration file for header update')
        parser.add_argument('--tilename', action='store', default=cls.DEFAULT_TILENAME,
                            help='Add (optional) TILENAME to SCI header')
        parser.add_argument('--tileid', action='store', type=int, default=cls.DEFAULT_TILEID,
                            help='Add (optional) TILENAME to SCI header')
        return

def get_safe_boolean(name,config,config_section):

    """ Get boolean first and if fail, get from config"""
    try:
        param = config.getboolean(config_section,name)
    except:
        param = config.get(config_section,name)
    return param

if __name__ == '__main__':
    CoaddZipperInterpNullWeight.main()
