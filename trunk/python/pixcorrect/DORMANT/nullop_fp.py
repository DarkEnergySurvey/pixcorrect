#!/usr/bin/env python
"""Do nothing to a DESImage
"""

# imports
from pixcorrect.PixCorrectDriver import PixCorrectFPStep

# constants

# Which section of the config file to read for this step
config_section = 'nullop_fp'

# exception classes
# interface functions
# classes

class NullOpFP(PixCorrectFPStep):
    description = "Do nothing"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """Execute a step that does nothing

        :Parameters:
            - `image`: the array of images to which to do nothing

        Applies the non-correction "in place"
        """
        return 0

nullop_fp = NullOpFP()

# internal functions & classes

if __name__ == '__main__':
    nullop_fp.main()
