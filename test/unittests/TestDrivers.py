"""Tests for command line drivers of pixcorrect tests
"""

import unittest
import logging
import logging.handlers
from unittest import TestCase, skip
from tempfile import mkdtemp
from os import path, environ, spawnv, P_WAIT, unlink, rmdir
from ConfigParser import SafeConfigParser
from contextlib import contextmanager
from shutil import rmtree

from TestPixCorrectIm import temp_pixcorrect_test_dir, ref_dir

LOG_FILENAME = 'TestDrivers.log'

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

    return logger

logger = prepare_logger()

class TestDrivers(TestCase):
    """Tests for command line drivers of pixcorrect tests

    These tests are not intended to test the functionality of the 
    code itself, just whether the code can be called from a 
    shell command.
    """
    sci_fname = path.join(ref_dir, 'scix.fits')
    bpm_fname = path.join(ref_dir, 'bpm.fits')

    def exec_runner(self, *args):
        prod_path = path.join(environ['PIXCORRECT_DIR'],'bin')
        cmd = path.join(environ['PIXCORRECT_DIR'],'bin',args[0])

        with temp_pixcorrect_test_dir() as temp_dir:
            out_fname = path.join(temp_dir, 'test_output.fits')
            args += ('-i', self.sci_fname, '-o', out_fname)
            logger.info('Running: ' + ' '.join(args))
            ret_code = spawnv(P_WAIT, cmd, args)

        self.assertEqual(ret_code, 0)

    def test_nullop(self):
        self.exec_runner('nullop')

    def test_mask_saturation(self):
        self.exec_runner('mask_saturation')

    def test_apply_bpm(self):
        self.exec_runner('apply_bpm','-b',self.bpm_fname)

    def test_override_bpm(self):
        self.exec_runner('override_bpm','-b',self.bpm_fname)

    def test_fix_cols(self):
        self.exec_runner('fix_cols','-b',self.bpm_fname)

    def test_pixcorrect_im(self):
        self.exec_runner('pixcorrect_im','--bpm',self.bpm_fname)

if __name__ == '__main__':
    unittest.main()
