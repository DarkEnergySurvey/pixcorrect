#!/usr/bin/env python

from pixcorrect.null_weights import null_weights
from pixcorrect.row_interp   import row_interp
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage
from despymisc.miscutils import elapsed_time
import time

import fitsio
import numpy as np

class RowInterpNullWeight(PixCorrectMultistep):

    config_section = "rowinterp_nullweight"
    description = 'Perform row_interp and null_weights in one step'
    step_name = config_section
    DEFAULT_ME_PREPARE = False

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

        # Check if we want special multi-epoch weighting
        me_prepare  = self.config.getboolean(self.config_section, 'me_prepare')
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
        return

if __name__ == '__main__':
    RowInterpNullWeight.main()
