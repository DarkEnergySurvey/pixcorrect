#!/usr/bin/env python
"""Do nothing to a DESImage
"""

# imports
from pixcorrect.PixCorrectDriver import PixCorrectImStep

# constants

# Which section of the config file to read for this step
config_section = 'nullop'

# exception classes
# interface functions
# classes

class NullOp(PixCorrectImStep):
    description = "Do nothing"
    step_name = config_section

    @classmethod
    def __call__(cls, image):
        """Execute a step that does nothing

        :Parameters:
            - `image`: the DESImage on which not to operate

        Applies the non-correction "in place"
        """
        return 0

nullop = NullOp()

# internal functions & classes

if __name__ == '__main__':
    nullop.main()
