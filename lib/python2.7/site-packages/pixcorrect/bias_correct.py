#!/usr/bin/env python
"""Example step bias in image correction"""

# imports
from functools import partial
import ctypes
import sys
from os import path

import numpy as np
import numpy.ma as ma

from pixcorrect import corr_util
from pixcorrect.harmonic_mean import harmonic_mean
from pixcorrect.corr_util import PrepMain, logger, do_once
from pixcorrect import proddir
from despyfits.DESImage import DESImage, DESDataImage

# constants

# Which section of the config file to read for this step
config_section = 'bias'

# exception classes
# interface functions

@do_once(0,'DESBIAS')
def bias_correct(image, config):
    """Apply operation bias correction to a DES image

    :Parameters:
        - `image`: the DESImage to be bias corrected
        - `config`: A python ConfigParser object

    Applies the correction "in place"
    """

    bias_fname = config.get(config_section, 'bias')
    bias_im = DESImage.load(bias_fname)
    image.header['BIASFIL'] = path.basename(bias_fname)
    image.data -= bias_im.data
    harmonic_mean(image.weight, bias_im.weight)

    logger.info('Completed bias')

# classes
# internal functions & classes

def bias_correct_main(config):
    """A driver for applying operation bias

    :Parameters:
        -`config`: a python ConfigParsers object

    @returns: 0 if successful
    """
    # Load files and images here if we don't want
    # bias to clean up the memory, because other drivers that
    # call bias might want it for future steps.
    in_fname = config.get(config_section, 'image_in_fname')
    im = DESImage.load(in_fname)

    bias_correct(im, config)

    out_fname = config.get(config_section, 'image_out_fname')
    im.save(out_fname, save_mask=False, save_weight=False)

    return 0

def add_bias_args(prep_main):
    prep_main.parser.add_argument('--image_in_fname', nargs=1, 
                                  default=None,
                                  help='input image file name')
    prep_main.parser.add_argument('--bias', nargs=1, 
                                  default=None, 
                                  help='bias filename')
    prep_main.parser.add_argument('--image_out_fname', nargs=1, 
                                  default=None,
                                  help='output image file name')
    prep_main.switch_opts['bias']=['image_in_fname', 
                                   'bias',
                                   'image_out_fname']

if __name__ == '__main__':
    prep_main = PrepMain(config_section, 'Apply a bias correction')
    add_bias_args(prep_main)
    config = prep_main()
    try:
        bias_correct_main(config)
        sys.exit(0)
    except:
        # If we want to set specific exit status values
        # based on what exceptions get thrown, do that
        # here
        raise
