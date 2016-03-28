#!/usr/bin/env python

from pixcorrect.null_weights import null_weights
from pixcorrect.row_interp   import row_interp
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage
from despymisc.miscutils import elapsed_time
from despyfits import updateWCS 
import time

import fitsio
import numpy as np

class CoaddRowInterpNullWeight(PixCorrectMultistep):

    config_section = "rowinterp_nullweight"
    description = 'Perform row_interp and null_weights in one step'
    step_name = config_section
    DEFAULT_ME_PREPARE = False
    DEFAULT_HEADFILE = False
    DEFAULT_HDUPCFG = False
    DEFAULT_TILENAME = False
    DEFAULT_TILEID = False

    # Fix the step_name for passing the command-line arguments to the classes
    null_weights.__class__.step_name = config_section
    row_interp.__class__.step_name   = config_section
    
    def __call__(self):
        """
        Run row_interp and null_weights in one step, we run the tasks
        by calling step_run in each class
        """
        t0 = time.time()
        # Get the science image
        input_image = self.config.get(self.config_section,'in')
        self.sci = DESImage.load(input_image)

        # Get verbose
        verbose = self.config.get(self.config_section,'verbose')

        # Check if we want special multi-epoch weighting
        me_prepare  = self.config.getboolean(self.config_section, 'me_prepare')

        # Get optional config file, first we try to get them as boolean, then as strings
        headfile = get_safe_boolean('headfile',self.config,self.config_section)
        hdupcfg  = get_safe_boolean('hdupcfg',self.config,self.config_section)

        self.update_sci_header(input_image)

        # Update the header if both headfile and hdupcfg are present
        if  headfile and hdupcfg:
            logger.info("Will update image header with scamp .head file %s" % headfile)
            self.sci = updateWCS.run_update(self.sci,headfile=headfile,hdupcfg=hdupcfg,verbose=verbose)
        
        if me_prepare:
            self.custom_weight(input_image)
        
        # Run null_weights
        t1 = time.time()
        logger.info("Running null_weights on: %s" % input_image)
        null_weights.step_run(self.sci,self.config)
        logger.info("Time NullWeights : %s" % elapsed_time(t1))

        # Run row_interp
        t2 = time.time()
        logger.info("Running row_interp on: %s" % input_image)
        row_interp.step_run(self.sci,self.config)
        logger.info("Time RowInterp : %s" % elapsed_time(t2))
        
        output_image = self.config.get(self.config_section, 'out')
        # Special write out
        if me_prepare:
            self.custom_write(output_image)
        else:
            self.sci.save(output_image)

        logger.info("Wrote new file: %s" % output_image)
        logger.info("Time Total: %s" % elapsed_time(t0))
        
        return 0

    def update_sci_header(cls,input_image):
        
        tilename = get_safe_boolean('tilename',cls.config,cls.config_section)
        tileid   = get_safe_boolean('tileid',cls.config,cls.config_section)
        if tilename:
            record={'name':'TILENAME', 'value':tilename, 'comment':'DES Tilename'}
            cls.sci.header.add_record(record)
        if tileid:
            record={'name':'TILEID', 'value':tileid, 'comment':'Tile ID for DES Tilename'}
            cls.sci.header.add_record(record)


    def custom_weight(cls,input_image):
        # Make custom weight, that will not zero STAR maskbit
        logger.info("Will perform special weighting for multi-epoch input on %s" % input_image)
        cls.weight_custom = np.copy(cls.sci.weight)
        null_mask = parse_badpix_mask(cls.config.get(cls.config_section, 'null_mask'))
        star_mask = parse_badpix_mask('STAR') # 32
        badamp_mask = parse_badpix_mask('BADAMP') 
        badamp    = np.array( cls.sci.mask & badamp_mask, dtype=bool)
        kill      = np.array( cls.sci.mask & null_mask, dtype=bool)
        stars     = np.array( cls.sci.mask & star_mask, dtype=bool)
        cls.weight_custom[kill]   = 0.0
        cls.weight_custom[stars]  = np.copy(cls.sci.weight[stars])
        cls.weight_custom[badamp] = 0.0

    def custom_write(cls,output_image):
        # Write out the image using fitsio, but skipping the mask as we won't need it.
        ofits = fitsio.FITS(output_image,'rw',clobber=True)
        ofits.write(cls.sci.data,  extname='SCI', header=cls.sci.header)
        ofits.write(cls.sci.weight,extname='WGT',header=cls.sci.weight_hdr)
        ofits.write(cls.weight_custom,extname='WGT_ME',header=cls.sci.weight_hdr)
        ofits.close()

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments for null_weights and row_interp
        """
        null_weights.add_step_args(parser)
        row_interp.add_step_args(parser)
        parser.add_argument('--me_prepare', action='store_true',default=cls.DEFAULT_ME_PREPARE,
                            help='Run custom weights for STAR and do not write MSK plane for multi-epoch (me)')
        parser.add_argument('--headfile', action='store', default=cls.DEFAULT_HEADFILE,
                            help='Headfile (containing most update information)')
        parser.add_argument('--hdupcfg', action='store', default=cls.DEFAULT_HDUPCFG,
                            help='Configuration file for header update')
        parser.add_argument('--tilename', action='store', default=cls.DEFAULT_TILENAME,
                            help='Add (optional) TILENAME to SCI header')
        parser.add_argument('--tileid', action='store', default=cls.DEFAULT_TILEID,
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
    RowInterpNullWeight.main()
