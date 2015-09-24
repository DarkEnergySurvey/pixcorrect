#!/usr/bin/env python

from pixcorrect.null_weights import null_weights
from pixcorrect.row_interp   import row_interp
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep

from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage, DESBPMImage

class RowInterpNullWeight(PixCorrectMultistep):

    config_section = "rinw"
    description = 'Run row_interp and null_weights in one step'
    step_name = config_section

    # Fixing the step name for pass the command-line arguments
    null_weights.__class__.step_name = config_section
    row_interp.__class__.step_name   = config_section
    
    def __call__(self):
        """
        Run row_interp and null_weights in one step'
        """

        # Get the science image
        input_image = self.config.get(self.config_section,'in')
        self.sci = DESImage.load(input_image)

        # Run the tasks by calling step_run in each class, passing the
        # triplet image and the config object
        logger.info("Running null_weights on: %s" % input_image)
        null_weights.step_run(self.sci,self.config)

        # Run row_interp
        logger.info("Running row_interp on: %s" % input_image)
        row_interp.step_run(self.sci,self.config)
        
        # Write out the image
        output_image = self.config.get(self.config_section, 'out')
        self.sci.save(output_image)
        logger.info("Wrote new file: %s" % output_image)

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
