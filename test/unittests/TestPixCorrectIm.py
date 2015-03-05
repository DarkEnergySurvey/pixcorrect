"""Tests for pixcorrect-im
"""

import unittest
import logging
import logging.handlers
from unittest import TestCase, skip, expectedFailure
from tempfile import mkdtemp
from os import path, environ
from ConfigParser import SafeConfigParser
from contextlib import contextmanager
from shutil import rmtree

import numpy as np

from pixcorrect.pixcorrect_im import PixCorrectIm
from despyfits.DESImage import DESImage
from despyfits import maskbits

LOG_FILENAME = 'TestPixCorrectIm.log'

def prepare_logger():
    logger = logging.getLogger(__name__)

    logger.setLevel(logging.DEBUG)

    old_log_exists = path.isfile(LOG_FILENAME)
    file_handler = logging.handlers.RotatingFileHandler(
        LOG_FILENAME, backupCount=99)
    file_handler.setLevel(logging.DEBUG)
    logfile_format = "%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s"
    logfile_formatter = logging.Formatter(logfile_format)
    file_handler.setFormatter(logfile_formatter)
    logger.addHandler(file_handler)
    if old_log_exists:
        logger.handlers[0].doRollover()

    # console_handler = logging.StreamHandler()
    # console_handler.setLevel(logging.DEBUG)
    # console_format = "%(message)s"
    # console_formatter = logging.Formatter(console_format)
    # console_handler.setFormatter(console_formatter)
    # logger.addHandler(console_handler)

    return logger

logger = prepare_logger()
ref_dir = environ['PIXCORRECTTESTDATA_DIR']

# Create a temp dir context manager to guarantee 
# the temp dir gets cleaned up, even if there are 
# exceptions thrown
#
# There is a standard function to do this in python 3.2,
# but I need to do it myself in 2.7
@contextmanager
def temp_pixcorrect_test_dir():
    temp_dir = mkdtemp(suffix='-pixcorrect-test')
    yield temp_dir
    rmtree(temp_dir, True)


# shorthand for adding a data file to the pixcorrect configuration
def add_data_config(config, keyword, fname, data_dir, section='pixcorrect_im'):
    if not section in config.sections():
        config.add_section(section)    
    fq_fname = path.join(data_dir, fname)
    config.set(section, keyword, fq_fname)

# shorthand for adding a data file in the reference directory 
# to the pixcorrect configuration
def add_ref_data_config(config, keyword, fname, section='pixcorrect_im'):
    add_data_config(config, keyword, fname, ref_dir, section)

class TestPixCorrectIm(TestCase):

    def new_config(self, tmp_dir):
        config = SafeConfigParser()
        add_ref_data_config(config, 'sci', 'scix.fits')
        add_data_config(config, 'out', 'corrected.fits', tmp_dir)
        config.set('pixcorrect_im', 'fixcol', 'False')
        return config
        
    def add_nullop_config(self, config):
        config.set('pixcorrect_im', 'nullop', 'True')

    def add_nullop_im_config(self, config):
        config.set('pixcorrect_im', 'nullop_im', 'True')

    def add_bias_config(self, config):
        add_ref_data_config(config, 'bias', 'biascor.fits')

    def add_flat_config(self, config):
        config.set('pixcorrect_im', 'flat', 'dflatcor.fits')

    def add_bpm_config(self, config):
        add_ref_data_config(config, 'bpm', 'bpm.fits')
        config.set('pixcorrect_im', 'mask_saturation', 'True')
        config.set('pixcorrect_im', 'fixcol', 'True')

    def add_override_bpm_config(self, config):
        self.add_bpm_config(config)
        add_ref_data_config(config, 'override_bpm', 'True')

    @expectedFailure
    def test_bias(self):
        with temp_pixcorrect_test_dir() as temp_dir:
            config = self.new_config(temp_dir)
            self.add_bpm_config(config)
            self.add_bias_config(config)
            pix_corrector = PixCorrectIm(config)
            logger.info('Doing bias correction')
            pix_corrector()

            test_im = DESImage.load( config.get('pixcorrect_im', 'out') )
            ref_im = DESImage.load( path.join(ref_dir, 'post_bias.fits') )
            in_im = DESImage.load( path.join(ref_dir, 'scix.fits') )
            im_cmp = ref_im.compare( test_im )
            logger.debug(str(im_cmp.header))
            im_cmp.log(logger, ref_im)
            self.assertTrue(im_cmp.match())

#    @skip("")
    def test_flat(self):
        with temp_pixcorrect_test_dir() as temp_dir:
            config = self.new_config(temp_dir)
            self.add_flat_config(config)
            pix_corrector = PixCorrectIm(config)
            logger.info('Doing flat correction')
            pix_corrector()

#    @skip("")
    @expectedFailure
    def test_bpm(self):
        with temp_pixcorrect_test_dir() as temp_dir:
            config = self.new_config(temp_dir)
            self.add_bpm_config(config)
            pix_corrector = PixCorrectIm(config)
            logger.info('Doing BPM correction')
            pix_corrector()

            test_im = DESImage.load( config.get('pixcorrect_im', 'out') )
            ref_im = DESImage.load( path.join(ref_dir, 'post_bpm.fits') )
            in_im = DESImage.load( path.join(ref_dir, 'scix.fits') )
            im_cmp = ref_im.compare( test_im )
            logger.debug(str(im_cmp.header))
            im_cmp.log(logger, ref_im)
            self.assertTrue(im_cmp.match())

#    @skip("")
#    @expectedFailure
    def test_nullop(self):
        with temp_pixcorrect_test_dir() as temp_dir:
            config = self.new_config(temp_dir)
            self.add_nullop_config(config)
            pix_corrector = PixCorrectIm(config)
            logger.debug('Doing nullop')
            pix_corrector()

    @expectedFailure
    def test_nullop_im(self):
        with temp_pixcorrect_test_dir() as temp_dir:
            config = self.new_config(temp_dir)
            self.add_nullop_im_config(config)
            pix_corrector = PixCorrectIm(config)
            logger.debug('Doing nullop_im')
            pix_corrector()

            test_im = DESImage.load( config.get('pixcorrect_im', 'out') )
            ref_im = DESImage.load( path.join(ref_dir, 'scix.fits') )
            im_cmp = ref_im.compare( test_im )
            logger.info("Nullop_Im results")
            logger.debug(str(im_cmp.header))
            im_cmp.log(logger, ref_im)
            self.assertTrue(im_cmp.match())

if __name__ == '__main__':
    unittest.main()
