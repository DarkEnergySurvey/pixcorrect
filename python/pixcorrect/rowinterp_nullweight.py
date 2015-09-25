#!/usr/bin/env python

from pixcorrect.null_weights import null_weights
from pixcorrect.row_interp   import row_interp
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage
from despymisc.miscutils import elapsed_time
import time

class RowInterpNullWeight(PixCorrectMultistep):

    config_section = "rowinterp_nullweight"
    description = 'Perform row_interp and null_weights in one step'
    step_name = config_section

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
        
        # Write out the image
        output_image = self.config.get(self.config_section, 'out')
        self.sci.save(output_image)
        logger.info("Wrote new file: %s" % output_image)
        logger.info("Time Total: %s" % elapsed_time(t0))
        
        return 0

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments for null_weights and row_interp
        """
        null_weights.add_step_args(parser)
        row_interp.add_step_args(parser)
        return

if __name__ == '__main__':
    RowInterpNullWeight.main()
