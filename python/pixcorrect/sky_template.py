#!/usr/bin/env python
"""
Use pre-established PCA coefficients to fit build full-res sky templates
Try a faster sigma-clipping approach that also uses mask images as a start
"""

from os import path
import numpy as np
import fitsio

from ConfigParser import SafeConfigParser, NoOptionError
from argparse import ArgumentParser
import time

from despyfits.DESImage import DESDataImage
from pixcorrect import proddir
from pixcorrect.corr_util import logger,items_must_match
from pixcorrect.PixCorrectDriver import PixCorrectDriver
from pixcorrect import skyinfo
from pixcorrect.clippedMean import clippedMean
from pixcorrect import decaminfo

# Which section of the config file to read for this step
config_section = 'skytemplate'

class SkyTemplate(PixCorrectDriver):
    description = "Create full-resolution sky templates based on previous PCA"
    step_name = config_section
    
    @classmethod
    def __call__(cls,
                 in_filename, out_filename,
                 ccdnum,
                 input_template=None,
                 input_list=None,
                 good_filename = None,
                 reject_rms = None,
                 mem_use = 8.,
                 bitmask = skyinfo.DEFAULT_SKYMASK):
        """
        Create full-resolution sky templates based on previous PCA.
        Does this pixel by pixel, via robust fitting of the data in the input
        full-res images to the PCA coefficients.  The full-res input images' filenames
        are determined from the EXPNUM _either_ by python formatting of the string
        given in -input_template _or_ by looking at the list of expnum, filename pairs
        in the file specified by -input-list.
        Output FITS image has an extension NGOOD giving number of
        images used in fit at each pixel.

        :Parameters:
            - `in_filename`: the file holding the PCA outputs on compressed sky
            - `out_filename`: filename for the output template
            - `ccdnum`: which CCD to produce templates for
            - `input_template`: string that can be formatted with the expnum to yield
                            filename of the DESImage holding the full-res data.
            - `input_list`: name of a file containing expnum, filename pairs, one pair per
                            line, separated by whitespace.
            - `good_filename`: Name of a FITS file in which to save number of images
                            contributing to each pixel's fit.  No output if None.  
            - `reject_rms`: Exclude exposures with fractional RMS residual sky above this.
                            If this is None, just uses the exposures that PCA used.
            - `mem_use:` Number of GB to target for memory usage (Default = 8)
            - `bitmask:` Applied to MASK extension of images for initial bad-pixel
                            exclusion.
        """
 
        logger.info('Starting sky template construction')

        # Need exactly one of these two arguments:
        if not ( (input_template is None) ^ (input_list is None)):
            logger.error('Need exactly one of input_template and input_list to be given')
            return 1
        
        # Acquire PCA information, including the table of info on input exposures
        pc = skyinfo.MiniskyPC.load(in_filename)
        pctab = skyinfo.MiniskyPC.get_exposures(in_filename)

        # Build a MiniDECam that has our choice of CCDs that we can use for indexing.
        mini = pc.get_pc(0)
        # Quit if we are requesting template for a CCD that was not compressed
        detpos = decaminfo.detpos_dict[ccdnum]
        try:
            mini.index_of(detpos,1,1)
        except skyinfo.SkyError:
            logger.error('Template requested for CCDNUM not included in PCA')
            return 1

        # Select exposures we'll use
        if reject_rms is None:
            # If no RMS threshold is specified, use the same exposures
            # that were kept during PCA of compressed skies
            use = np.array(pctab['USE'])
        else:
            # Choose our own threshold
            use = pctab['RMS'] < reject_rms

        # Get filenames for the full-res images from list:
        if input_list is not None:
            filenames = {}
            flist = np.loadtxt(input_list,dtype=str)
            for expnum,filename in flist:
                filenames[int(expnum)] = filename
            del flist
            # Now warn if we are missing expnums and remove from usable exposure list
            for i in range(len(use)):
                if use[i] and not int(pctab['EXPNUM'][i]) in filenames.keys():
                    use[i] = False
                    logger.warning('No input filename given for expnum ' + str(expnum))
                    
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

        # The linear solutions will require this:
        ainv = np.linalg.inv( np.dot(V,V.T))
        
        nexp = V.shape[1]
        npc = pc.U.shape[1]
        ySize = decaminfo.shape[0]
        xSize = decaminfo.shape[1]

        # Create the output array
        out = np.zeros( (npc, ySize, xSize), dtype=np.float32)

        # And an array to hold the number of exposures used at each pixel:
        if good_filename is not None:
            ngood = np.zeros( (ySize, xSize), dtype=np.int16)

        # Decide how many rows of blocks we'll read from files at a time
        bytes_per_row = 4 * xSize * pc.blocksize * nimg
        xBlocks = xSize / pc.blocksize
        yBlocks = min( int(np.floor( 0.8* mem_use * (2**30) / bytes_per_row)),
                       ySize / pc.blocksize)

        if yBlocks < 1:
            logger.warning('Proceeding even though mem_use is not enough to store 1 row of blocks')
            yBlocks = 1
            
        d = {'ccd':ccdnum}   # A dictionary used to assign names to files
        hdr = {}      # A dictionary of information to go into output image header
        # A mask of zero is equivalent to no masking:
        if bitmask==0:
            bitmask = None
            
        nonConvergentBlocks = 0  # Keep count of blocks where clipping does not converge.
        
        # Collect input data in chunks of yBlocks rows of blocks, then process one block at a time.
        for yStart in range(0,ySize,yBlocks*pc.blocksize):
            # Acquire the pixel data into a 3d array
            yStop = min(ySize,yStart+yBlocks*pc.blocksize)
            logger.info('Working on rows {:d} -- {:d}'.format(yStart,yStop))
            data = np.zeros( (nimg, yStop-yStart, xSize), dtype=np.float32)
            # Mask image:
            mask = np.zeros( (nimg, yStop-yStart, xSize), dtype=bool)

            for i,expnum in enumerate(expnums):
                d['expnum']=expnum
                if input_template is None:
                    # Get the filename from the input list
                    filename = filenames[expnum]
                else:
                    # Get the filename from formatting the template
                    filename = input_template.format(**d)
                logger.debug('Getting pixels from ' + filename)
                with fitsio.FITS(filename) as fits:
                    data[i,:,:] = fits['SCI'][yStart:yStop, :xSize]
                    if bitmask is None:
                        mask[i,:,:] = True
                    else:
                        m = np.array(fits['MSK'][yStart:yStop, :xSize],dtype=np.int16)
                        mask[i,:,:] = (m & bitmask)==0
                        del m
                    if yStart==0:
                        # First time through the images we will be collecting/checking
                        # header information from the contributing images
                        hdrin = fits['SCI'].read_header()
                        usehdr = {}
                        if 'BAND' in hdrin.keys():
                            usehdr['BAND'] = hdrin['BAND']
                        elif 'FILTER' in hdrin.keys():
                            usehdr['BAND'] = decaminfo.get_band(hdrin['FILTER'])
                        else:
                            logger.error('No BAND or FILTER in ' + filename)
                            return 1
                        if 'NITE' in hdrin.keys():
                            usehdr['NITE'] = hdrin['NITE']
                        elif 'DATE-OBS' in hdrin.keys():
                            usehdr['NITE'] = decaminfo.get_nite(hdrin['DATE-OBS'])
                        else:
                            logger.error('No NITE or DATE-OBS in ' + filename)
                            return 1
                        if 'FLATFIL' in hdrin.keys():
                            usehdr['FLATFIL'] = hdrin['FLATFIL']
                        else:
                            logger.error('No FLATFIL in ' + filename)
                            return 1
                        if 'CCDNUM' in hdrin.keys():
                            usehdr['CCDNUM'] = hdrin['CCDNUM']
                        else:
                            logger.error('No CCDNUM in ' + filename)
                            return 1
                        if len(hdr)==0:
                            # First exposure will establish values for the output
                            hdr['BAND'] = usehdr['BAND']
                            hdr['MINNITE'] = usehdr['NITE']
                            hdr['MAXNITE'] = usehdr['NITE']
                            hdr['CCDNUM'] = usehdr['CCDNUM']
                            if hdr['CCDNUM']!=ccdnum:
                                logger.error('Wrong ccdnum {:d} in {:s}'.format(
                                    ccdnum,filename))
                            hdr['FLATFIL'] = usehdr['FLATFIL']
                        else:
                            # Check that this exposure matches the others
                            try:
                                items_must_match(hdr, usehdr, 'BAND','CCDNUM','FLATFIL')
                            except:
                                return 1
                            hdr['MINNITE'] = min(hdr['MINNITE'],usehdr['NITE'])
                            hdr['MAXNITE'] = max(hdr['MAXNITE'],usehdr['NITE'])
                            
            data /= norms[:,np.newaxis,np.newaxis]  # Apply norms to be near unity
                    
            # Now cycle through all blocks
            for jb in range((yStop-yStart)/pc.blocksize):
                for ib in range(xSize/pc.blocksize):
                    logger.debug('Fitting for block ({:d},{:d})'.format(jb+yStart/pc.blocksize,ib))
                    if ccdnum==decaminfo.ccdnums['S7'] and \
                       pc.halfS7 and \
                       ib >= xSize/pc.blocksize/2:
                        # If we are looking at the bad amp of S7, we'll just
                        # store the median of the normalized images in PC0.
                        # The other PC's stay at zero.
                        out[0,
                            yStart+jb*pc.blocksize:yStart+(jb+1)*pc.blocksize,
                            ib*pc.blocksize:(ib+1)*pc.blocksize] = \
                          np.median(data[:,
                                         jb*pc.blocksize:(jb+1)*pc.blocksize,
                                         ib*pc.blocksize:(ib+1)*pc.blocksize],
                                    axis=0)
                        continue
            

                    # Use PCA of this block as starting guess at solution
                    index = mini.index_of(detpos,
                                          yStart/pc.blocksize + jb,
                                          ib)
                    guess = np.array(pc.U[index,:])

                    # Extract the data for this block into (nexp,npix) array
                    block = np.array(data[:,
                                        jb*pc.blocksize:(jb+1)*pc.blocksize,
                                        ib*pc.blocksize:(ib+1)*pc.blocksize])
                    block.resize(nexp, pc.blocksize*pc.blocksize)

                    bmask = np.array(mask[:,
                                          jb*pc.blocksize:(jb+1)*pc.blocksize,
                                          ib*pc.blocksize:(ib+1)*pc.blocksize])
                    bmask.resize(nexp, pc.blocksize*pc.blocksize)

                    # We'll scale the guess in each pixel by the typical ratio
                    # of this pixel's data to the PCA model for the block, and
                    # also estimate noise as dispersion about this guess
                    model = np.dot(guess, V)
                    ratio = block / model[:,np.newaxis]
                    scale, var, n = clippedMean(ratio,4,axis=0)
                    clip = 3. * np.sqrt(var.data)*scale.data

                    # First guess at solution is the outer product of superblock PCA
                    # with the scaling per pixel
                    soln = guess[:,np.newaxis]*scale.data
                    del scale, var, ratio, n

                    # Linear solution with clipping iteration
                    MAX_ITERATIONS = 20
                    TOLERANCE = 0.0001
                    for i in range(MAX_ITERATIONS):
                        model = np.dot(V.T,soln)
                        # Residuals from model are used to clip
                        resid = block - model
                        # Find clipped points and masked ones
                        good = np.logical_and(resid < clip, resid > -clip)
                        good = np.logical_and(good, bmask)
                        # Set residual to zero at bad pixels
                        resid[~good] = 0.
                        # Get shift in linear solution from residuals:
                        dsoln = np.dot(ainv, np.dot(V,resid))
                        soln += dsoln
                        # Calculate largest change in model as convergence criterion
                        shift = np.max(np.abs(np.dot(V.T,dsoln)))
                        logger.debug('Iteration {:d}, model shift {:f}'.format(i,shift))
                        if shift < TOLERANCE:
                            break
                        if i==MAX_ITERATIONS-1:
                            nonConvergentBlocks = nonConvergentBlocks + 1
                                           
                    # Save results into big matrices
                    soln.resize(npc,pc.blocksize,pc.blocksize)
                    out[:,
                        yStart+jb*pc.blocksize:yStart+(jb+1)*pc.blocksize,
                        ib*pc.blocksize:(ib+1)*pc.blocksize] = soln
                    if good_filename is not None:
                        # Gin up a masked array because it allows counting along an axis
                        nblock = np.ma.count_masked(\
                            np.ma.masked_array(np.zeros_like(good),good),axis=0)
                        nblock.resize(pc.blocksize,pc.blocksize)
                        ngood[yStart+jb*pc.blocksize:yStart+(jb+1)*pc.blocksize,
                              ib*pc.blocksize:(ib+1)*pc.blocksize] = nblock
                        del nblock
                    del resid, model, good, dsoln, block
            del data

        if nonConvergentBlocks > 0:
            logger.warning('Clipping did not converge for {:d} blocks out of {:d}'.format\
                           (nonConvergentBlocks,xBlocks*(ySize/pc.blocksize)))

        # Add a history line about creation here
        hdr['HISTORY'] = time.asctime(time.localtime()) + \
            ' Build sky template from PCA file {:s}'.format(path.basename(in_filename))
        # Save the template into the outfile
        spc = skyinfo.SkyPC(out,detpos,header=hdr)
        spc.save(out_filename)
        del out
        
        # Save the number of good sky pixels in another extension
        if good_filename is not None:
            gimg = DESDataImage(ngood, header={'DETPOS':detpos,
                                               'CCDNUM':ccdnum})
            logger.debug('Writing ngood to ' + good_filename)
            gimg.save(good_filename)
            del gimg, ngood
            
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
        mem_use = config.getfloat(cls.step_name,'mem_use')

        if config.has_option(cls.step_name,'input_template'):
            input_template = config.get(cls.step_name,'input_template')
        else:
            input_template = None

        if config.has_option(cls.step_name,'input_list'):
            input_list = config.get(cls.step_name,'input_list')
        else:
            input_list = None
            
        if config.has_option(cls.step_name,'reject_rms'):
            reject_rms = config.getfloat(cls.step_name,'reject_rms')
        else:
            reject_rms = None

        if config.has_option(cls.step_name,'good_filename'):
            good_filename = config.get(cls.step_name,'good_filename')
        else:
            good_filename = None

        ret_code = cls.__call__(in_filename=infile,
                                out_filename=out_filename,
                                ccdnum=ccdnum,
                                input_template=input_template,
                                input_list=input_list,
                                reject_rms=reject_rms,
                                mem_use=mem_use,
                                good_filename=good_filename,
                                bitmask = skyinfo.DEFAULT_SKYMASK)
        return ret_code

    @classmethod
    def add_step_args(cls,parser):

        parser.add_argument('-i','--infile',type=str,
                            help='File with PCA information (from sky_pca)')
        parser.add_argument('-o','--outfilename',type=str,
                            help='Name for output FITS template file')
        parser.add_argument('-c','--ccdnum',type=int,
                            help='CCDNUM of device for which to build templates')
        parser.add_argument('--input_template',type=str,
                            help='String which yields filenames of individual FITS images when formatted,'
                            ' e.g. D{expnum:08d}_{ccd:02d}_fp.fits')
        parser.add_argument('--input_list',type=str,
                            help='File holding pairs of expnum, filename on each line giving input images')
        parser.add_argument('--reject_rms', type=float,
                            help='Reject exposures with RMS resids from PCA fit above this')
        parser.add_argument('--good_filename', type=str,
                            help='FITS file to hold counts of valid exposures per pixel')
        parser.add_argument('--mem_use', type=float, default=8.,
                            help='Number of GB of memory usage to target')
        return


sky_template = SkyTemplate()

# internal functions & classes

if __name__ == '__main__':
    sky_template.main()
