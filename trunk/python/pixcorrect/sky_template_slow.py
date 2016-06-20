#!/usr/bin/env python
"""
Use pre-established PCA coefficients to fit build full-res sky templates
"""

from os import path
import numpy as np
import fitsio

from ConfigParser import SafeConfigParser, NoOptionError
from argparse import ArgumentParser

from pixcorrect import proddir
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectImStep
from pixcorrect import skyinfo
from pixcorrect.clippedMean import clippedMean
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'skytemplate'

class SkyTemplate(PixCorrectImStep):
    description = "Create full-resolution sky templates based on previous PCA"
    step_name = config_section
    
    @classmethod
    def __call__(cls, in_filename, out_filename, ccdnum, img_template, reject_rms, mem_use):
        """
        Create full-resolution sky templates based on previous PCA.
        Does this pixel by pixel, via robust fitting of the data in the input
        full-res images to the PCA coefficients.

        :Parameters:
            - `in_filename`: the file holding the PCA outputs on compressed sky
            - `out_filename`: filename for the output template
            - `ccdnum`: which CCD to produce templates for
            - `img_template`: string that can be formatted with the expnum to yield
                              filename of the DESImage holding the full-res data.
            - `reject_rms`: Exclude exposures with fractional RMS residual sky above this.
                            If this is None, just uses the exposures that PCA used.
            - `mem_use:` Number of GB to target for memory usage
        """
 
        logger.info('Starting sky template construction')

        # Acquire PCA information, including the table of info on input exposures
        pc = skyinfo.MiniskyPC.load(in_filename)
        # ??? Should make this table part of the MiniskyPC class:
        pctab = fitsio.read(in_filename,ext='EXPOSURES')

        # Build a MiniDECam that has our choice of CCDs that we can use for indexing.
        mini = pc.get_pc(0)
        # Quit if we are requesting template for a CCD that was not compressed
        detpos = decaminfo.detpos_dict[ccdnum]
        try:
            mini.index_of(detpos,1,1)
        except SkyError:
            logger.error('Template requested for CCDNUM not included in PCA')
            return(1)

        # Select exposures we'll use
        if reject_rms is None:
            # If no RMS threshold is specified, use the same exposures
            # that were kept during PCA of compressed skies
            use = np.array(pctab['USE'])
        else:
            # Choose our own threshold
            use = pctab['RMS'] < reject_rms
        nimg = np.count_nonzero(use)

        expnums = []
        vv = []
        for i in range(len(use)):
            if use[i]:
                vv.append(pctab['COEFFS'][i])
                expnums.append(pctab['EXPNUM'][i])
        V = np.vstack(vv)
        del vv
        
        # We'll re-normalize each exposure, and its coefficients, by V[0]
        norms = np.array(V[:,0])
        V = V.T/norms  # V is now of shape (npc,nimg)

        npc = pc.U.shape[1]
        ySize = decaminfo.shape[0]
        xSize = decaminfo.shape[1]

        # Create the output array
        out = np.zeros( (npc, ySize, xSize), dtype=np.float32)

        # Only fill half of it for the bad amp:
        if ccdnum==decaminfo.ccdnums['S7'] and pc.halfS7:
            xSize = xSize/2
            
        # Decide how many rows of blocks we'll read from files at a time
        bytes_per_row = 4 * xSize * pc.blocksize * nimg
        xBlocks = xSize / pc.blocksize
        yBlocks = min( int(np.floor( mem_use * (2**30) / bytes_per_row)),
                       ySize / pc.blocksize)

        if yBlocks < 1:
            logger.warning('Proceeding even though mem_use is not enough to store 1 row of blocks')
            yBlocks = 1
            
        d = {'ccd':ccdnum}

        # Collect input data in chunks of yBlocks rows of blocks, then process one block at a time.
        for yStart in range(0,ySize,yBlocks*pc.blocksize):
            # Acquire the pixel data into a 3d array
            yStop = min(ySize,yStart+yBlocks*pc.blocksize)
            logger.info('Working on rows {:d} -- {:d}'.format(yStart,yStop))
            data = np.zeros( (nimg, yStop-yStart, xSize), dtype=np.float32)

            for i,expnum in enumerate(expnums):
                d['expnum']=expnum
                filename = img_template.format(**d)
                logger.debug('Getting pixels from ' + filename)
                with fitsio.FITS(filename) as fits:
                    data[i,:,:] = fits['SCI'][yStart:yStop, :xSize]
            data /= norms[:,np.newaxis,np.newaxis]  # Apply norms to be near zero
                    
            # Now cycle through all blocks
            for jb in range((yStop-yStart)/pc.blocksize):
                for ib in range(xSize/pc.blocksize):
                    logger.debug('Fitting for block ({:d},{:d})'.format(jb+yStart/pc.blocksize,ib))
                    # Use PCA of this block as starting guess at solution
                    index = mini.index_of(detpos,
                                          yStart/pc.blocksize + jb,
                                          ib)
                    guess = np.array(pc.U[index,:])

                    # We'll scale the guess in each pixel by the typical ratio
                    # of this pixel's data to the PCA model for the block:
                    model = np.dot(guess, V)
                    ratio = data[:,
                                 jb*pc.blocksize:(jb+1)*pc.blocksize,
                                 ib*pc.blocksize:(ib+1)*pc.blocksize] / model[:,np.newaxis,np.newaxis]
                    scale, var, n = clippedMean(ratio,4,axis=0)
                    logger.debug('Var, scale, ratio shapes: ' + str(var.shape) + \
                                 ' ' + str(scale.shape) + ' ' + str(ratio.shape))

                    del ratio, n
                    # Solve each pixel in the block:
                    for jp in range(pc.blocksize):
                        for ip in range(pc.blocksize):
                            cost = skyinfo.ClippedCost(3*np.sqrt(var[jp,ip]))
                            # Execute and save the fit
                            out[:, yStart+jb*pc.blocksize+jp, ib*pc.blocksize+ip] = \
                              skyinfo.linearFit(data[:, jb*pc.blocksize+jp, ib*pc.blocksize+ip],
                                                V,
                                                guess*scale[jp,ip],
                                                cost)

            del data

        # Save the template into the outfile
        spc = skyinfo.SkyPC(out,detpos)
        spc.save(out_filename, clobber=True)
        
        logger.debug('Finished sky template')
        ret_code=0
        return ret_code

    @classmethod
    def run(cls, config):
        """Customized execution for sky template.  

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """

        infile = config.get(cls.step_name, 'infile')
        out_filename = config.get(cls.step_name, 'outfilename')
        ccdnum = config.getint(cls.step_name, 'ccdnum')
        reject_rms = config.getfloat(cls.step_name, 'reject_rms')
        mem_use = config.getfloat(cls.step_name,'mem_use')
        image_template = config.get(cls.step_name,'image_template')
        if config.has_option(cls.step_name,'reject_rms'):
            reject_rms = config.getfloat(cls.step_name,'reject_rms')
        else:
            reject_rms = None

        ret_code = cls.__call__(in_filename=infile,
                                out_filename=out_filename,
                                ccdnum=ccdnum,
                                img_template=image_template,
                                reject_rms=reject_rms,
                                mem_use=mem_use)
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

        parser.add_argument('-i','--infile',type=str,
                            help='File with PCA information (from sky_pca)')
        parser.add_argument('-o','--outfilename',type=str,
                            help='Name for output FITS template file')
        parser.add_argument('-c','--ccdnum',type=int,
                            help='CCDNUM of device for which to build templates')
        parser.add_argument('--image_template',type=str,default='D{expnum:08d}_{ccd:02d}_fp.fits',
                            help='String which yields filenames of individual FITS images when formatted')
        parser.add_argument('--reject_rms', type=float,
                            help='Reject exposures with RMS resids from PCA fit above this')
        parser.add_argument('--mem_use', type=float, default=8.,
                            help='Number of GB of memory usage to target')
        return parser


sky_template = SkyTemplate()

# internal functions & classes

if __name__ == '__main__':
    sky_template.main()
