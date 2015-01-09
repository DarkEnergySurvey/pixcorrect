"""Tests for foo
"""

from unittest import TestCase, skip
from tempfile import mkdtemp
from os import path
from ConfigParser import SafeConfigParser
from contextlib import contextmanager
from shutil import rmtree

import numpy as np
import pyfits

import pixcorrect.foo as foo
from pixcorrect.imtypes import ImageWrongHeader

# Create a temp dir context manager to guarantee 
# the temp dir gets cleaned up, even if there are 
# exceptions thrown
#
# There is a standard function to do this in python 3.2,
# but I need to do it myself in 2.7
@contextmanager
def temp_foo_test_dir():
    temp_dir = mkdtemp(suffix='TestFoo')
    yield temp_dir
    rmtree(temp_dir, True)


class TestFoo(TestCase):
    imsize = (2048, 4096)
    random_seed = 6563
    coeff = 42.02

    def make_config(self, im2_fname):
        config = SafeConfigParser()
        config.add_section('foo')
        config.set('foo', 'im2_fname', im2_fname)
        config.set('foo', 'im2_hdu', '0')
        config.set('foo', 'foo_coeff', "%f" % self.coeff)
        return config

    def make_hdu1(self):
        data1 = 2*np.random.random(self.imsize)-1
        hdu1 = pyfits.PrimaryHDU(data1)
        hdu1.header.set('FORALL', True)
        hdu1.header.set('TREE', 'oak')
        return hdu1

    def make_hdu2(self, temp_dir):
        im2_fname = path.join(temp_dir, 'im2.fits')
        data2 = 2*np.random.random(self.imsize)-1
        hdu2 = pyfits.PrimaryHDU(data2)
        hdu2.header['FORFOO'] = True
        hdu2.header['BOGUS'] = 'nope'
        hdulist = pyfits.HDUList([hdu2])
        hdulist.writeto(im2_fname)
        return data2, im2_fname
        
    def test_foo(self):
        hdu1 = self.make_hdu1()
        ar1_orig = hdu1.data.copy()
        with temp_foo_test_dir() as d:
            ar2, im2_fname = self.make_hdu2(d)
            config = self.make_config(im2_fname)
            ret = foo.foo(hdu1, config)

        # Test the results
        expected = ar1_orig*self.coeff + ar2
        difference = hdu1.data-expected
        maxdiff = np.amax(difference)
        self.assertEqual(maxdiff, 0)
        mindiff = np.amin(difference)
        self.assertEqual(mindiff, 0)

    def test_foo_badhdu1(self):
        data1 = 2*np.random.random(self.imsize)-1
        hdu1 = pyfits.PrimaryHDU(data1)
        hdu1.header.set('TREE', 'oak')

        ar1_orig = hdu1.data.copy()
        with temp_foo_test_dir() as d:
            ar2, im2_fname = self.make_hdu2(d)
            config = self.make_config(im2_fname)
            self.assertRaises(ImageWrongHeader,
                              foo.foo, hdu1, config)
