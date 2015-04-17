#!/usr/bin/env python
"""
Perform robust PCA on ensemble of mini-sky images and save results to file.
"""

from os import path
import numpy as np
from ConfigParser import SafeConfigParser, NoOptionError
from argparse import ArgumentParser

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo
from pixcorrect.clippedMean import clippedMean

# Which section of the config file to read for this step
config_section = 'skypca'

def rank1(m):
    """
    Create an initial rank-1 model for data m by taking
    medians along rows and columns iteratively.
    Returns
    med, amps
    which are the median image (normed to unit median) and amplitudes
    per exposure.  Outer product of these is the model (med[:,np.newaxis]*amps).
    """
    # Iterate medians a couple of times:
    med = np.ones(m.shape[0],dtype=float)
    for i in range(2):
        amps = np.median(m/med[:,np.newaxis],axis=0)
        med = np.median(m/amps,axis=1)
    # Set median image to have unit median:
    norm = np.median(med)
    med /= norm
    amps *= norm
    return med,amps

def pca(m,nkeep=20):
    """
    Perform PCA on data matrix m.  Return
    u*s,s,v
    for the top nkeep PC's.
    """
    R = np.dot(m.T, m)
    s,v = np.linalg.eigh(R)
    # Sort eigenvalues and retain just top nkeep
    idx = np.argsort(s)[::-1][:nkeep]
    v = v[:,idx]
    s = s[idx]
    return np.dot(v.T, m.T).T,np.sqrt(s),v

def clip(data,model,nsigma=4):
    diff = data - model
    avg,var,n = clippedMean(diff,nsigma,maxSample=1000000)
    print 'Mean = ',avg,'+-',np.sqrt(var)
    diff = np.abs(diff-avg) > nsigma*np.sqrt(var)
    out = np.where(diff, model+avg, data)
    return out, np.count_nonzero(diff)

def process(m, npc):
    """
    Sequence of operations to perform a robust PCA
    """
    # Construct rank-1 approximation to data and scale it out
    # so all entries are near unity
    med, amps = rank1(m)
    data = m / med[:,np.newaxis]
    data /= amps
        
    # The initial model is all 1's
    model = np.ones_like(data)
    # Now iterate clip, PCA process
    npc_clip = max(npc, 10)   # Number of PCs used during clipping
    n_clip_iterations = 3   # Rounds of clipping of outliers
    n_clip_sigma = 4        # Number of sigma for clipping during iteration
    n_clip_sigma_final = 3  # Final clipping pass
    for i in range(n_clip_iterations):
        work,n = clip(data,model,n_clip_sigma)
        us,s,v = pca(work,npc)
        model = np.dot(us[:,:npc],v[:,:npc].T)

    # And a final pass with tighter clip
    work,n = clip(data,model,n_clip_sigma_final)
    us,s,v = pca(work,npc)
    
    # Restore the rank-1 scalings
    U = us[:,:npc] * med[:,np.newaxis]
    V = v[:,:npc] * amps[:,np.newaxis]
    # Now rescale the us to have unit variance along columns
    norms = np.sqrt(np.sum(U*U,axis=0)/U.shape[0])
    U /= norms
    V *= norms
    return U, V, s

class SkyPCA(PixCorrectImStep):
    description = "Perform robust PCA of a collection of MiniDecam images"
    step_name = config_section
    
    @classmethod
    def __call__(cls, in_filenames, out_filename, npc, reject_rms):
        """
        Perform robust PCA of a collection of MiniDecam images.

        :Parameters:
            - `in_filenames`: list of filenames of exposure mini-sky image
            - `out_filename`: filename for the output PCA
            - `npc`: Number of principal components to retain
            - `reject_rms`: Exclude exposures with fractional RMS above this
        """
 
        logger.info('Collecting images for PCA')
        mm = []  # will hold the data vectors for each exposure
        data_length = None
        blocksize = None
        for f in in_filenames:
            mini = skyinfo.MiniDecam.load(f)
            v = mini.vector()
            if data_length is None:
                data_length = len(v)
            elif len(v) != data_length:
                logger.error('Mismatched sky data vector length in file ' + f)
                raise SkyError('Mismatched sky data vector length in file ' + f)
            # Also check for agreement of the mini-sky setup
            if blocksize is None:
                blocksize = mini.blocksize
                mask_value = mini.mask_value
                invalid = mini.invalid
                print invalid, "read in"
                halfS7 = mini.halfS7
            else:
                if mini.blocksize != blocksize \
                  or mini.invalid!=invalid \
                  or mini.halfS7!=halfS7:
                  logger.error('Mismatched minisky configuration in file ' + f)
                  raise SkyError('Mismatched minisky configuration in file ' + f)
            mm.append(np.array(v))
        m = np.vstack(mm).transpose()
        del mm

        logger.info("Start first PCA cycle")
        U, S, v = process(m, npc)
        
        pc = skyinfo.MiniskyPC(U,
                               blocksize=blocksize,
                               mask_value=mask_value,
                               invalid = invalid,
                               halfS7 = halfS7)

        # Refit each exposure to this PCA,
        nexp = len(in_filenames)
        V = np.zeros((nexp, U.shape[1]),dtype=float)
        rms = np.zeros(nexp,dtype=float)
        frac = np.zeros(nexp,dtype=float)
        for i in range(nexp):
            mini.fill_from(m[:,i])
            pc.fit(mini,clip_sigma=3.)
            rms[i] = mini.rms
        use = rms <= reject_rms
        logger.info('Retained {:d} out of {:d} exposures'.format(np.count_nonzero(use),
                                                                 nexp))

        # New PCA excluding outliers
        logger.info("Start second PCA cycle")
        U, S, v = process(m, npc)
        
        pc = skyinfo.MiniskyPC(U,
                               blocksize=blocksize,
                               mask_value=mask_value,
                               invalid = invalid,
                               halfS7 = halfS7)

        pc.save(out_filename,clobber=True)
        
        # Recollect statistics and save
        logger.info("Collecting statistics")
        for i in range(nexp):
            mini.fill_from(m[:,i])
            pc.fit(mini,clip_sigma=3.)
            V[i,:] = mini.coeffs
            rms[i] = mini.rms
            frac[i] = mini.frac

        # ??? Write V and a results table using fitsio

        # Create a one-line binary fits table to hold the coefficients
        logger.debug('Finished PCA')
        ret_code=0
        return ret_code

    @classmethod
    def run(cls, config):
        """Customized execution for sky combination.  Note there is NO input image nor output

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        inlist = open(config.get(cls.step_name, 'inlist'))
        in_filenames = []
        for line in inlist:
            in_filenames.append(line.strip())
        out_filename = config.get(cls.step_name, 'outfilename')
        npc = config.getint(cls.step_name, 'npc')
        reject_rms = config.getfloat(cls.step_name, 'reject_rms')

        ret_code = cls.__call__(in_filenames, out_filename, npc, reject_rms)
        return ret_code

    @classmethod
    def parser(cls):
        """Generate a parser
        """
        default_config = path.join(proddir, 'etc', cls.step_name+'.config')
        default_out_config = path.join(cls.step_name+'-as_run'+'.config')

        # Argument parser
        parser = ArgumentParser(description=cls.description)
        parser.add_argument("config", default=default_config, nargs="?",
                            help="Configuration file filename")
        parser.add_argument('-s', '--saveconfig', 
                                 default=default_out_config,
                                 help="output config file")
        parser.add_argument('-l', '--log', 
                                 default=cls.step_name+".log", 
                                 help="the name of the logfile")
        parser.add_argument('-v', '--verbose', action="count", 
                                 help="be verbose")

        parser.add_argument('-i','--inlist',type=str,
                            help='File holding names of all input minisky files')
        parser.add_argument('-o','--outfilename',type=str,
                            help='Filename for PCA results/resids')
        parser.add_argument('-n','--npc',type=int,default=4,
                            help='Number of principal components to retain')
        parser.add_argument('--reject_rms', type=float, default=0.005,
                            help='Reject exposures with RMS resids from PCA fit above this')
        return parser


sky_pca = SkyPCA()

# internal functions & classes

if __name__ == '__main__':
    sky_pca.main()
