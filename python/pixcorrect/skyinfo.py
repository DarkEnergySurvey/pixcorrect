"""
Classes and defaults for sky subtraction
"""

import copy
import numpy as np
import fitsio
from despyfits.DESImage import data_dtype
from despyfits import maskbits
from pixcorrect import decaminfo
from pixcorrect import clippedMean

DEFAULT_BLOCKSIZE = 128  # Decimation factor for sky images
DEFAULT_BITMASK = \
    maskbits.BADPIX_BPM | \
    maskbits.BADPIX_SATURATE | \
    maskbits.BADPIX_EDGE  # bitmask for pixels to ignore in sky calculations
DEFAULT_IGNORE = 'N30,S30,S7'  # Chips to leave out of all sky calculations
DEFAULT_MASK_VALUE = -1.  # Value assigned to unspecified pixels in compressed sky image
DEFAULT_CCDNUMS = '1-62'  # CCDNUMs to use
DEFAULT_MINISKY_FILES = 'minisky{:>02d}.fits'  # Template for mini-sky file names to assemble.
DEFAULT_CLIP_SIGMA = 3.  # Rejection threshold for robust fitting of sky to PCs.
DEFAULT_SKY_WEIGHT = True # Whether to build weight image from fitted sky image

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

    `header` property is a FITS header that is used to store other properties, or
    any supplied keyword/value information.

    `coeffs` property is array of coefficients of sky templates for fit to this. Stored
    in the header as 'SKYPC00','SKYPC01',...

    `rms` property is fractional RMS residuals to fit.  Uses keyword 'SKYRMS'

    `frac` property is fraction of outliers in sky fit.  Uses keyword 'SKYFRAC'
    """

    def __init__(self,
                 blocksize=DEFAULT_BLOCKSIZE,
                 mask_value=DEFAULT_MASK_VALUE,
                 invalid=('N30','S30','S7'),
                 header = None):
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
        self.header = fitsio.FITSHDR(header)

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

    @property
    def coeffs(self):
        """
        Array of coefficients for sky templates.  Stored in header as 'SKYPCnn' entries
        """
        c = []
        MAX_PC = 20
        for ipc in range(MAX_PC):
            kw = 'SKYPC{:>02d}'.format(ipc)
            if kw not in self.header:
                break
            c.append(self.header[kw])
        return np.array(c, dtype=float)

    @coeffs.setter
    def coeffs(self,c):
        for ipc,val in enumerate(c):
            kw = 'SKYPC{:>02d}'.format(ipc)
            self.header[kw] = float(val)
        return
        
    @property
    def rms(self):
        """
        RMS fractional deviation of sky from model, after sigma-clipping
        """
        return self.header['SKYRMS']

    @rms.setter
    def rms(self,rms):
        self.header['SKYRMS'] = rms
        return

    @property
    def frac(self):
        """
        Fraction of valid sky pixels rejected by sigma clipping when calculating RMS
        """
        return self.header['SKYFRAC']

    @frac.setter
    def frac(self,frac):
        self.header['SKYFRAC'] = frac
        return
        
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

    def save(self,filename):
        """
        Save the mini-image to primary extension of a FITS file.
        """
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
        d,hdr = fitsio.read(filename,header=True)
        blocksize = hdr['BLOCKSIZ']
        mask_value = hdr['MASKVAL']
        invalid = hdr['BADDET'].split(',')
        out = cls(blocksize, mask_value, invalid, header=hdr)
        out.data = d
        return out

class MiniskyPC(object):
    """
    Class containing principle components of compressed sky images.
    """
    def __init__(self, U):
        """
        :Parameters:

        `U`: 2d numpy array with shape (number of sky superpixels : number of pcs)
        """
        self.U = np.array(U)
        return

    @classmethod
    def load(cls,filename,ext='U'):
        """
        Retrieve PC's from a FITS file, from primary or specified extension
        """
        U = fitsio.read(filename,ext=ext)
        return cls(U)

    def save(self,filename,ext='U', clobber=True, header=None):
        """
        Save the PCs into a FITS file of the given name, in specified extension
        Clobber=True by default; otherwise will always append a new extension to the file
        """
        fitsio.write(filename, self.U, extname=ext,clobber=clobber, header=header)
        return

    def fit(self, mini, clip_sigma=3.):
        """
        Determine coefficients of PCs that best fit the given MiniDecam.  Iteratively fit
        using a robust fitter which downweights deviant pixels, presumably due to object
        in the sky.

        Upon return, `mini` image is replaced by fractional residuals to the fit, and
        the statistics of the fit are stored in it.

        :Parameters:

            - `mini`:  a MiniDecam image to be fit.
            - `clipSigma`: Number of sigma to mark as outliers in robust fitting
        """

        # Extract the data vector from the mini-image
        y = mini.vector()
        # Restrict the data and the templates to the superpixels where we have valid data
        use = y != mini.mask_value
        y = y[use]
        x = self.U[use,:]
        
        # Create first guess as simply the 0th template, scaled to match median
        aStart = np.zeros(x.shape[1],dtype=float)
        aStart[0] = np.median(y / x[:,0])
        # Determine a sigma for the residuals and build a cost function, 4-sigma clipping
        avg,var,n = clippedMean(y-x[:,0]*aStart[0],4)
        cost = ClippedCost(4*np.sqrt(var))
        # Initial fit
        a = linearFit(y, x.T, aStart, cost)
        # Repeat fit with updated variance estimate and clipping threshold
        avg,var,n = clippedMean(y-np.dot(x,a),4)
        cost = ClippedCost(clip_sigma*np.sqrt(var))
        a = linearFit(y, x.T, a, cost)

        # Get statistics of fractional residuals to this fit
        avg,var,n = clippedMean(y/np.dot(x,a)-1.,clip_sigma)
        v = mini.vector()
        v[use] = y / np.dot(x,a) - 1.
        mini.fill_from(v)
        mini.coeffs = a
        mini.rms = np.sqrt(var)
        mini.frac = 1.-float(n)/len(data)
        return

class SkyPC(object):
    """
    Full-resolution sky principal components (templates)
    """
    # ??? Add ability to use templates that have been subsampled to save I/O time.

    def __init__(self, d, detpos):
        """
        A sky pc is a 3d array with index 0 enumerating 2d principal components of sky.
        """
        self.d = np.array(d)
        self.detpos = detpos
        return

    @classmethod
    def load(cls, filename, ext=None):
        """
        Get a sky pc from the 3d array stored in named extension of FITS file
        """
        d,h = fitsio.read(filename, ext=ext, header=True)
        if len(d.shape) != 3:
            raise SkyError("SkyTemplates.load did not find 3d array in " + filename)
        detpos = h['DETPOS']
        return cls(d,detpos)

    def save(self, filename, d, ext=None, clobber=True, header=None):
        """
        Save a sky pc as a FITS file under given extension name (or primary).
        If clobber=False, it is appended as a new extension.
        """
        if h is None:
            h = {}
        else:
            h = fitsio.FITSHDR(header)
        h['DETPOS'] = self.detpos
        fitsio.write(filename, self.d, extname=ext,clobber=clobber, header=h)
        return
    def sky(self, coeffs):
        """
        Return a 2d array constructed from applying the specified coefficients
        """
        if len(coeffs.shape)!=1 or self.d.shape[0] != len(coeffs):
            raise SkyError("Wrong number of coefficients for SkyTemplates.sky: " +
                           str(coeffs.shape))
        return np.sum(coeffs[:,np.newaxis,np.newaxis]*self.d, axis=0)

def linearFit(y, x, aStart, cost, dump=False):
    """
    Find the value of vector a that minimizes the sum of cost
    function applied to the elements of vector
    y - ax
    y and x are given data vector and coefficient matrix.
    Cost must be a class that returns the total value, first, and second derivs
    of the cost (i.e. negative log likelihood)
    given the residual vector as input.

    This will be solved iteratively for a, beginning with aStart. Convergence judged
    on change in the cost function.

    Returns the final a.  Returns all zeros if does not converge.
    """
    a = np.array(aStart)
    oldCost = None
    iterations = 0
    MAX_ITERATIONS = 10
    COST_TOLERANCE = 0.01
    while (iterations < MAX_ITERATIONS):
        dy = y - np.dot(a, x)
        totcost, d1, d2 = cost(dy)
        if dump:
            print iterations,':',totcost, a
        if oldCost!=None:
            if oldCost - totcost < COST_TOLERANCE*oldCost:
                # done!
                return a
        oldCost = totcost
            
        # do an update
        beta = np.dot(x, d1)
        alpha = np.dot(x*d2, x.T)
        a += np.linalg.solve(alpha, beta)
    print "Too many iterations in linearFit"
    return np.zeros_like(aStart)

class ClippedCost:
    """
    A cost function that implements least-squares fitting in linearFit with
    clipping beyond the prescribed value of limit.
    """
    def __init__(self, limit):
        self.limit = limit
        return
    def __call__(self, data):
        use = np.abs(data) < self.limit
        cost = np.sum((data*data)[use]) / np.count_nonzero(use)
        d1 = np.where(use,data,0.)
        d2 = np.where(use, 1., 0.)
        return cost,d1,d2

