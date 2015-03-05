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
    fpsci_fname = path.join(ref_dir, 'DECam_00394250.fits.fz')
    fpout_base_fname = 'test_output_%%d.fits'

    sci_fname = path.join(ref_dir, 'scix.fits')
    out_base_fname = 'test_output.fits'
    bpm_fname = path.join(ref_dir, 'bpm.fits')
    bias_fname = path.join(ref_dir, 'biascor.fits')
    flat_fname = 'dflatcor.fits'

    def exec_runner(self, in_fname, out_base_fname, *args):
        prod_path = path.join(environ['PIXCORRECT_DIR'],'bin')
        cmd = path.join(environ['PIXCORRECT_DIR'],'bin',args[0])

        with temp_pixcorrect_test_dir() as temp_dir:
            out_fname = path.join(temp_dir, out_base_fname)
            args += ('-i', in_fname, '-o', out_fname)
            logger.info('Running: ' + ' '.join(args))
            ret_code = spawnv(P_WAIT, cmd, args)

        self.assertEqual(ret_code, 0)

    def im_exec_runner(self, *args):
        self.exec_runner(self.sci_fname, self.out_base_fname, *args)

    def fp_exec_runner(self, *args):
        self.exec_runner(self.fpsci_fname, self.fpout_base_fname, *args)

#    @skip('')
    def test_nullop(self):
        prod_path = path.join(environ['PIXCORRECT_DIR'],'bin')
        cmd = path.join(environ['PIXCORRECT_DIR'],'bin','nullop')

        args = (cmd, )
        with temp_pixcorrect_test_dir() as temp_dir:
            logger.info('Running: ' + ' '.join(args))
            ret_code = spawnv(P_WAIT, cmd, args)

        self.assertEqual(ret_code, 0)

#    @skip('')
    def test_nullop_im(self):
        self.im_exec_runner('nullop_im')

#    @skip('')    
    def test_bias_correct(self):
        self.im_exec_runner('bias_correct','--bias', self.bias_fname)

#    @skip('')    
    def test_mask_saturation(self):
        self.im_exec_runner('mask_saturation')

#    @skip('')    
    def test_apply_bpm(self):
        self.im_exec_runner('apply_bpm','-b', self.bpm_fname)

#    @skip('')    
    def test_override_bpm(self):
        self.im_exec_runner('override_bpm','-b', self.bpm_fname)

#    @skip('')    
    def test_fix_cols(self):
        self.im_exec_runner('fix_cols','-b', self.bpm_fname)

#    @skip('')    
    def test_flat_correct(self):
        self.im_exec_runner('flat_correct','--flat', self.flat_fname)

#    @skip('')    
    def test_pixcorrect_im(self):
        self.im_exec_runner('pixcorrect_im','--bpm',self.bpm_fname)

#    @skip('')    
    def test_pixcorrect_fp(self):
        self.fp_exec_runner('pixcorrect_fp','--nullop_fp')

if __name__ == '__main__':
    unittest.main()
