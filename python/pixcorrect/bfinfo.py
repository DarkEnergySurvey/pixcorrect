#!/usr/bin/env python
"""
Data structure used to hold the brighter/fatter corrections
"""

import numpy as np
import fitsio
from despyfits.DESImage import data_dtype
from despyfits import maskbits

DEFAULT_BFMASK = \
    maskbits.BADPIX_BPM | \
    maskbits.BADPIX_SATURATE | \
    maskbits.BADPIX_TRAIL | \
    maskbits.BADPIX_EDGEBLEED | \
    maskbits.BADPIX_BADAMP # bitmask for pixels to ignore in shift calculations

class BFError(Exception):
    """
    Exception class for problems in brighter-fatter correction
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class BFKernel(object):
    """
    Class that holds the coefficients for brighter/fatter relation
    """
    
    def __init__(self, bffile, detpos):
        """
        Read brighter/fatter coefficients from a saved FITS file.  The upper-right
        quadrant information for right (top) edge shifts is stored in FITS image
        extensions named e.g. "N12R", "N12T"

        :Parameters:
           - `bffile`: path to the FITS file holding coefficient arrays
           - `detpos`: detector for which to seek information.
        """
        ff = fitsio.FITS(bffile)
        if detpos + 'R' not in ff:
            raise BFError('Kernel table "' + detpos + 'R' + '" not in BF file ' + bffile)
        self.aR = ff[detpos+'R'].read()

        if detpos + 'T' not in ff:
            raise BFError('Kernel table "' + detpos + 'T' + '" not in BF file ' + bffile)
        self.aT = ff[detpos+'T'].read()

        return

    def kernelR(self, shape):
        """
        Return an array of specified shape that has the full 2d kernel filled in.
        1 electron located at [0,0] causes the specified fractional shift of the
        *right* edge of each pixel in the array.  Pixels to left/below of the charge
        wrap around the array.
        """
        if shape[0] < 2*self.aR.shape[0]-1 or shape[1] < 2*self.aR.shape[1]+1:
            raise BFError('Requested shape for kernelR is too small to hold the kernel')
        kernel = np.zeros(shape, dtype=data_dtype)
        ny, nx = self.aR.shape
        #r_kernel will contain the fraction of pixel to shift out of right edge
        kernel[:ny,:nx] = self.aR
        # Replicate the other quadrants accordingly
        kernel[-1:-ny:-1,:nx] = self.aR[1:,:]
        kernel[:ny, -1:-nx-1:-1] = -self.aR
        kernel[-1:-ny:-1,-1:-nx-1:-1] = -self.aR[1:,:]
        return kernel

    def kernelT(self, shape):
        """
        Return an array of specified shape that has the full 2d kernel filled in.
        """
        if shape[0] < 2*self.aT.shape[0]+1 or shape[1] < 2*self.aT.shape[1]-1:
            raise BFError('Requested shape for kernelT is too small to hold the kernel')
        kernel = np.zeros(shape, dtype=data_dtype)
        ny, nx = self.aT.shape
        kernel[:ny,:nx] = self.aT
        # Replicate the other quadrants accordingly
        kernel[:ny, -1:-nx:-1,] = self.aT[:,1:]
        kernel[-1:-ny-1:-1, :nx] = -self.aT
        kernel[-1:-ny-1:-1, -1:-nx:-1] = -self.aT[:,1:]
        return kernel
    
    @classmethod
    def from_gruen(cls, gruenfile, outfile):
        """
        Translate from Daniel's format into FITS format.
        """
        fin = open(gruenfile)
        fout = fitsio.FITS(outfile,'rw',clobber=True)
        for line in fin:
            fields = line.split()
            detpos = fields.pop(0).strip()
            if fields.pop(0) != 'M':
                raise BFError('Did not get "M" in Gruen-format BF file')
            nx = int(fields.pop(0))
            ny = int(fields.pop(0))
            aR = np.zeros( (ny,nx), dtype=data_dtype)
            for ix in range(nx):
                for iy in range(ny):
                    aR[iy,ix] = np.float32(fields.pop(0))
            if fields.pop(0) != 'M':
                raise BFError('Did not get "M" in Gruen-format BF file')
            nx = int(fields.pop(0))
            ny = int(fields.pop(0))
            aT = np.zeros( (ny,nx), dtype=data_dtype)
            for ix in range(nx):
                for iy in range(ny):
                    aT[iy,ix] = np.float32(fields.pop(0))
            if len(fields)!=0:
                raise BFError('Too much info in bf coefficient line for detpos ' + detpos)
            ## Axis swap...needed according to Daniel:
            aR,aT = aT.transpose(), aR.transpose()
            fout.write(aR, extname=detpos+'R')
            fout.write(aT, extname=detpos+'T')
        return
    
