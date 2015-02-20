"""
Classes and defaults for sky subtraction
"""

import numpy as np
import fitsio
from despyfits.DESImage import data_dtype
from despyfits import maskbits
from pixcorrect import decaminfo

DEFAULT_BLOCKSIZE = 128  # Decimation factor for sky images
DEFAULT_BITMASK = \
    maskbits.BADPIX_BPM | \
    maskbits.BADPIX_SATURATE | \
    maskbits.BADPIX_EDGE  # bitmask for pixels to ignore in sky calculations
DEFAULT_IGNORE = 'N30,S30,S7'  # Chips to leave out of all sky calculations
DEFAULT_MASK_VALUE = -1.  # Value assigned to unspecified pixels in compressed sky image
DEFAULT_CCDNUMS = '1-62'  # CCDNUMs to use
DEFAULT_MINISKY_FILES = 'minisky{:>02d}.fits'  # Template for mini-sky file names to assemble.

class SkyError(Exception):
    """
    Error class for problems in sky-subtraction methods
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def parse_ranges(ccdlist):
    """
    Take a string of format like "1,3,5-9,23-42" and return a list of
    integers that are specified.
    Raise ValueError if the unexpected is encountered.  Whitespace is ok.
    """
    s = set()
    for r1 in ccdlist.split(','):
        r2 = r1.split('-')
        if len(r2)==1:
            set.add(int(r2[0]))
        elif len(r2)==2:
            for j in range(int(r2[0]),int(r2[1])+1):
                s.add(j)
        else:
            raise ValueError('Bad integer range expression in parse_ranges: ' + r1)
    return sorted(s)

            

class MiniDecam(object):
    """
    Class holding a decimated image of the full DECam science array.
    Each pixel in this image represents a (blocksize x blocksize) region in
    the original array.

    The 2d array is accessed via 'data' property.

    `mask` property is a 2d boolean array which is True in pixels that fall on
    the CCDs that are in use.

    """

    def __init__(self,
                 blocksize=DEFAULT_BLOCKSIZE,
                 mask_value=DEFAULT_MASK_VALUE,
                 invalid=('N30','S30','S7')):
        """
        blocksize is the size of superpixels used in decimation of the data.  Must be
        a power of 2 to equally divide the DECam CCDs

        maskvalue will be inserted into all pixels that are not filled with valid data.

        invalid is a list of DETPOS values for CCDs that will not be used.
        """
        self.blocksize = blocksize
        self.mask_value = mask_value
        self.invalid = set()
        if invalid is not None:
            for detpos in invalid:
                self.invalid.add(detpos.strip())

        self._chip = [decaminfo.shape[0] / blocksize,
                       decaminfo.shape[1] / blocksize]  # The shape of a decimated CCD
        if decaminfo.shape[0]%self._chip[0] != 0 or \
           decaminfo.shape[1]%self._chip[1] != 0:
            # Raise exception if image is not multiple of blocksize.
            raise SkyError('MiniImage blocksize ' + str(blocksize) +
                            ' does not evenly divide images')

        # Get the bounding box for CCDs in use
        self.xmin = 0
        self.ymin = 0
        x0 = None
        x1 = None
        y0 = None
        y1 = None
        for detpos in decaminfo.ccdnums.keys():
            if detpos in self.invalid:
                continue
            y,x = self._corner_of(detpos)
            if x0 is None:
                x0 = x
                x1 = x + self._chip[1]
                y0 = y
                y1 = y + self._chip[0]
            x0 = min(x0,x)
            x1 = max(x1, x+self._chip[1])
            y0 = min(y0, y)
            y1 = max(y1, y+self._chip[0])
        self.xmin = x0
        self.ymin = y0

        # Create the data and mask images
        self.data = np.ones( (y1-y0,x1-x0), dtype=data_dtype) * self.mask_value
        self.mask = np.zeros(self.data.shape, dtype=bool)

        # Mark all useful regions in mask:
        for detpos in decaminfo.ccdnums.keys():
            if detpos in self.invalid:
                continue
            y,x = self._corner_of(detpos)
            self.mask[y:y+self._chip[0], x:x+self._chip[1]] = True
                
        return

    def _corner_of(self,detpos):
        """
        Return 2d coordinates of the (0,0) pixel of this detector.
        """
        if detpos in self.invalid or detpos not in decaminfo.ccdnums.keys():
            raise SkyError('Invalid detpos in MiniDecam: ' + detpos)
            
        x = (decaminfo.ccdCorners[detpos][0]-1-self.xmin)/self.blocksize
        y = (decaminfo.ccdCorners[detpos][1]-1-self.ymin)/self.blocksize
        return y,x
            
    def fill(self,data, detpos):
        """
        Fill the portion of the mini-image corresponding to detpos with the
        array given by data.
        """
        y,x = self._corner_of(detpos)
        self.data[ y:y+self._chip[0], x:x+self._chip[1] ] = data
        return

    def vector(self):
        """
        Return a 1d array, a flattened version of the array that contains only the
        pixels that are on valid CCDs.
        """
        return self.data[self.mask]

    def fill_from(self,vectorIn):
        """
        Set the data array equal to the values in the flattened array vectorIn
        """
        self.data[self.mask] = vectorIn
        return

    def index_of(self,detpos, j, i):
        """
        Return the index in the flattened vector() that would contain pixel (j,i)
        in *compressed* version of the chosen CCD.  (j,i) are in the numpy convention.
        Raises an exception for an invalid detpos.
        """
        tmp = np.zeros_like(self.mask)
        y,x = self._corner_of(detpos)
        tmp[y+j, x+i] = True
        return np.where(tmp[self.mask])[0][0]

    def save(self,filename, header=None):
        """
        Save the mini-image to primary extension of a FITS file.
        header gives desired header information to attach.
        """
        if header is None:
            header = {}
        header['BLOCKSIZ'] = self.blocksize
        header['MASKVAL'] = self.mask_value
        baddet = ''
        for detpos in self.invalid:
            if len(baddet)>0:
                baddet = baddet + ','
            baddet = baddet + detpos
        header['BADDET'] = baddet
        fitsio.write(filename, self.data, header=header, clobber=True)
        return

    @classmethod
    def load(cls, filename):
        """
        Extracts mini-image and header from primary extension of the given FITS file.
        Returns image,header
        """
        d,h = fitsio.read(filename,header=True)
        blocksize = h['BLOCKSIZ']
        mask_value = h['MASKVAL']
        invalid = h['BADDET'].split(',')
        out = cls(blocksize, mask_value, invalid)
        out.data = d
        return out,h


