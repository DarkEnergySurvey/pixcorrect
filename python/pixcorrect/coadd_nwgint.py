#!/usr/bin/env python3

import time
import os
import re
import sys

import matplotlib.path
import fitsio
import numpy as np

from pixcorrect.null_weights import null_weights
from pixcorrect.row_zipper   import row_zipper
from pixcorrect.corr_util import logger
from pixcorrect.PixCorrectDriver import PixCorrectMultistep
from pixcorrect.clip_mask_utils import polygon_to_pix

from despyastro.CCD_corners import update_DESDM_corners
from despyastro import wcsutil, astrometry

import despyfits
from despyfits.maskbits import parse_badpix_mask
from despyfits.DESImage import DESImage, update_hdr_compression, insert_eupspipe
from despyfits import updateWCS

from despymisc.miscutils import elapsed_time


class CoaddZipperInterpNullWeight(PixCorrectMultistep):

    """Run custom weights for STAR and do not write MSK plane for multi-epoch (me)'"""

    config_section = "coadd_nwgit"
    description = 'Perform zipper interpolation along rows and null_weights in one step'
    step_name = config_section
    DEFAULT_CUSTOM_WEIGHT = False
    DEFAULT_HEADFILE = False
    DEFAULT_HDUPCFG = False
    DEFAULT_TILENAME = False
    DEFAULT_TILEID = False
    DEFAULT_ME_WGT_KEEPMASK = False
    DEFAULT_NULL_MASK_SCI = '0'
#   region_mask additions
    DEFAULT_SKIP_IMGCHK = False

    # Fix the step_name for passing the command-line arguments to the classes
    null_weights.__class__.step_name = config_section
    row_zipper.__class__.step_name = config_section

    def __call__(self):
        """
        Run row_zipper and null_weights in one step, we run the tasks
        by calling step_run in each class
        """
        t0 = time.time()

        # Check if we want special multi-epoch weighting, and which bits we want to 'save'
        me_wgt_keepmask = get_safe_boolean('me_wgt_keepmask', self.config, self.config_section)

        # Get verbose
        try:
            verbose = self.config.get(self.config_section, 'verbose')
        except:
            verbose = False

        # Get the science image
        input_image = self.config.get(self.config_section, 'in')
        self.sci = DESImage.load(input_image)

        # In case a streak table is provided -- we proceed with the extra STREAK masking
        streak_file = self.config.get(self.config_section, 'streak_file')
        if os.path.exists(streak_file):
            add_width = self.config.getfloat(self.config_section, 'add_width')
            add_length = self.config.getfloat(self.config_section, 'add_length')
            max_extrapolate = self.config.getfloat(self.config_section, 'max_extrapolate')
            self.streakMask(streak_file,
                            addWidth=add_width,
                            addLength=add_length,
                            maxExtrapolate=max_extrapolate)

        # Add TILENAME and TILEID to sci header (optional) if required
        self.update_sci_header(input_image)

        # Update the header wcs if both headfile and hdupcfg are present (optional)
        self.update_wcs_header(input_image, verbose=verbose)

        # Case region mask is provided -- we proceed with the extra region masking
        region_mask = self.config.get(self.config_section, 'region_mask')
        if (os.path.exists(region_mask)):
            reg_flag_val = self.config.get(self.config_section, 'reg_flag_val')
            reg_skip_imgchk = get_safe_boolean('region_skip_imgchk', self.config, self.config_section)
            self.regionMask(region_mask,reg_flag_val=reg_flag_val,skip_expccd_chk=reg_skip_imgchk)

        # Check if want to create the custon weight for SWArp/SExtractor combination
        if me_wgt_keepmask:
            self.custom_weight(input_image)

        # Run null_weights
        t1 = time.time()
        logger.info("Running null_weights on: %s", input_image)
        null_weights.step_run(self.sci, self.config)
        logger.info("Time NullWeights : %s", elapsed_time(t1))

        # Run row_zipper
        t2 = time.time()
        logger.info("Running row_zipper on: %s", input_image)
        row_zipper.step_run(self.sci, self.config)
        logger.info("Time ZipperInterp : %s", elapsed_time(t2))

        # Null the sci image only if null_mask_sci !=0
        self.null_sci(input_image)

        output_image = self.config.get(self.config_section, 'out')
        # Special write out
        if me_wgt_keepmask:
            self.custom_write(output_image)
        else:
            self.sci.save(output_image)

        logger.info("Wrote new file: %s", output_image)
        logger.info("Time Total: %s", elapsed_time(t0))

        return 0

    def update_wcs_header(self, input_image, verbose=False):

        # Get optional config file, first we try to get them as boolean, then as strings
        headfile = get_safe_boolean('headfile', self.config, self.config_section)
        hdupcfg = get_safe_boolean('hdupcfg', self.config, self.config_section)

        # Update the header if both headfile and hdupcfg are present
        if  headfile and hdupcfg:
            logger.info("Will update image header with scamp .head file %s", headfile)
            self.sci = updateWCS.run_update(self.sci, headfile=headfile, hdupcfg=hdupcfg, verbose=verbose)

    def update_sci_header(self, input_image):
        tilename = get_safe_boolean('tilename', self.config, self.config_section)
        tileid = get_safe_boolean('tileid', self.config, self.config_section)
        if tilename:
            record = {'name': 'TILENAME', 'value': tilename, 'comment': 'DES Tilename'}
            self.sci.header.add_record(record)
        if tileid:
            record = {'name': 'TILEID', 'value': int(tileid), 'comment': 'Tile ID for DES Tilename'}
            self.sci.header.add_record(record)


    def null_sci(self, input_image):

        null_mask_sci = parse_badpix_mask(self.config.get(self.config_section, 'null_mask_sci'))
        if null_mask_sci != 0:
            logger.info('Nulling science image from null_mask_bits')
            kill = np.array(self.sci.mask & null_mask_sci, dtype=bool)
            self.sci.data[kill] = 0.0
        else:
            logger.info('Science image was not null')

    def custom_weight(self, input_image):
        # Make custom weight, that will not zero STAR maskbit
        logger.info("Will perform special weighting for multi-epoch input on %s", input_image)
        # Make a copy of the original untouched weight
        self.sci.weight_custom = np.copy(self.sci.weight)
        null_mask = self.config.get(self.config_section, 'null_mask')
        me_wgt_keepmask = self.config.get(self.config_section, 'me_wgt_keepmask')

        # Make python lists of the coma-separated input lists
        null_list = null_mask.split(',')
        keep_list = me_wgt_keepmask.split(',')

        # Special case we care:
        # . we are nulling the TRAIL but want keep where STAR
        if 'TRAIL' in null_list and 'STAR' in keep_list and 'TRAIL' not in keep_list:
            # Remove STAR from the list
            if 'STAR' in null_list:
                null_list.remove('STAR')
            null_mask_bits = parse_badpix_mask(','.join(null_list))
            # Null each plane at a time. First the TRAILS and replace with STAR
            kill = np.array(self.sci.mask & parse_badpix_mask('TRAIL'), dtype=bool)
            stars = np.array(self.sci.mask & parse_badpix_mask('STAR'), dtype=bool)
            self.sci.weight_custom[kill] = 0.0
            self.sci.weight_custom[stars] = np.copy(self.sci.weight[stars])
            # Loop over the bitplanes, but skipping TRAIL, which we already did
            null_list.remove('TRAIL')
            for bitplane in null_list:
                kill = np.array(self.sci.mask & parse_badpix_mask(bitplane), dtype=bool)
                self.sci.weight_custom[kill] = 0.0
        # We  remove tham from the null_list
        else:
            for bitplane in me_wgt_keepmask.split(','):
                if bitplane in null_list:
                    null_list.remove(bitplane)
            null_mask_bits = parse_badpix_mask(','.join(null_list))
            kill = np.array(self.sci.mask & null_mask_bits, dtype=bool)
            self.sci.weight_custom[kill] = 0.0

    def custom_write(self, output_image):
        # Write out the image using fitsio, but skipping the mask as we won't need it.
        ofits = fitsio.FITS(output_image, 'rw', clobber=True)

        # Here we mimick the steps followed by DESImage.save()
        # SCI
        logger.info("Creating SCI HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        self.sci.header = update_hdr_compression(self.sci.header, 'SCI')
        logger.info("Calculating CCD corners/center/extern keywords for SCI HDU ")
        self.sci.header = update_DESDM_corners(self.sci.header, get_extent=True, verb=False)
        if despyfits.DESImage.pipekeys_write:
            logger.info("Inserting EUPS PIPEPROD and PIPEVER to SCI HDU")
            self.sci.header = insert_eupspipe(self.sci.header)
        ofits.write(self.sci.data, extname='SCI', header=self.sci.header)
        # WGT
        logger.info("Creating WGT HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        self.sci.weight_hdr = update_hdr_compression(self.sci.weight_hdr, 'WGT')
        ofits.write(self.sci.weight, extname='WGT', header=self.sci.weight_hdr)
        # WGT_ME
        # For  WGT_ME we do not need to update the FZ keywords, as we use the same hdr as WGT
        logger.info("Creating WGT_ME HDU")
        ofits.write(self.sci.weight_custom, extname='WGT_ME', header=self.sci.weight_hdr)
        # MSK
        logger.info("Creating MSK HDU and relevant FZ*/DES_EXT/EXTNAME keywords")
        self.sci.mask_hdr = update_hdr_compression(self.sci.mask_hdr, 'MSK')
        ofits.write(self.sci.mask, extname='MSK', header=self.sci.mask_hdr)
        ofits.close()


    def streakMask(self, streak_file, addWidth=0., addLength=100., maxExtrapolate=0):
        '''
        Produce a list of pixels in the image that should be masked for
        streaks in the input table.  streaktab is the output table of new
        streaks to add image is a FITS HDU, with header and image data
        addWidth is additional number of pixels to add to half-width
        addLength is length added to each end of streak (pixels)

        Returns:
        ypix, xpix: 1d arrays with indices of affected pixels
        nStreaks: number of new streaks masked
        '''

        # Read the streaks table first
        try:
            tab = fitsio.FITS(streak_file)
            streaktab = tab[1].read()
        except:
            logger.error('Could not read streak file {:s}'.format(streak_file))
            sys.exit(1)

        image_header = self.sci.header
        image_data = self.sci.data
        # Pixscale in degrees
        pixscale = astrometry.get_pixelscale(image_header, units='arcsec') / 3600.
        shape = image_data.shape

        # # Due to a bug in fitsio 1.0.0rc1+0, we need to clean up the
        # # header before feeding it to wcsutil and remove the 'None' and other problematic items
        # for k in image_header:
        #     # Try to access the item, if failed we hace to remove it
        #     try:
        #         item = image_header[k]
        #     except:
        #         logger.info("Removing keyword: {:s} from header".format(k))
        #         image_header.delete(k)

        w = wcsutil.WCS(image_header)

        # WE NEED TO UPDATE THIS WHEN THE TABLE IS PER EXPNUM
        use = np.logical_and(streaktab['expnum'] == image_header['EXPNUM'],
                             streaktab['ccdnum'] == image_header['CCDNUM'])
        logger.info('{:d} streaks found to mask'.format(np.count_nonzero(use)))

        nStreaks = 0
        inside = None


        for row in streaktab[use]:
            if maxExtrapolate > 0:
                if row['extrapolated'] and row['nearest'] > maxExtrapolate:
                    logger.info('Skipping extrapolated streak')
                    continue
            width = row['width']
            ra = np.array((row['ra1'], row['ra2']))
            dec = np.array((row['dec1'], row['dec2']))
            x, y = w.sky2image(ra, dec)

            x1, x2, y1, y2 = x[0], x[1], y[0], y[1]

            # Slope of the line, cos/sin form
            mx = (x2 - x1) / np.hypot(x2 - x1, y2 -y1)
            my = (y2 - y1) / np.hypot(x2 - x1, y2 -y1)

            #displacement for width of streak:
            wx = width / pixscale + addWidth
            wy = wx * mx
            wx = wx * -my

            # grow length
            x1 -= addLength * mx
            x2 += addLength * mx
            y1 -= addLength * my
            y2 += addLength * my

            # From Alex's immask routine: mark interior pixels
            vertices = [(x1 + wx, y1 + wy), (x2 + wx, y2 + wy), (x2 - wx, y2 - wy), (x1 - wx, y1 - wy)]
            vertices.append(vertices[0])  # Close the path

            if inside is None:
                # Set up coordinate arrays
                yy, xx = np.indices(shape)
                points = np.vstack((xx.flatten(), yy.flatten())).T
                path = matplotlib.path.Path(vertices)
                inside = path.contains_points(points)
            else:
                # use logical_and for additional streaks
                path = matplotlib.path.Path(vertices)
                inside = np.logical_or(inside, path.contains_points(points))

            nStreaks = nStreaks + 1

        logger.info('Masked {:d} new streaks'.format(nStreaks))

        # Make the list of masked pixels
        if inside is None:
            ymask, xmask = np.array(0, dtype=int), np.array(0, dtype=int)
        else:
            ymask, xmask = np.nonzero(inside.reshape(shape))

        logger.info('Setting bits in MSK image for STREAK: {:d}'.format(parse_badpix_mask('STREAK')))
        self.sci.mask[ymask, xmask] |= parse_badpix_mask('STREAK')



    def regionMask(self, region_file, reg_flag_val="CRAY", skip_expccd_chk=False):
        '''Read a region file and mask pixels based on sets of points and polygons
        The current assumption for polygons is that they are convex (i.e. a convex hull)
        like those generated by the clip_mask utility (implemented to analyze outlier pixels
        inforation from SWarp COMBINE_TYPE=CLIPPED).
        
        Inputs: 
            region file: (basically a DS9 region file expressing regions in SKY coordinates)
            reg_flag_val: (default = "CRAY")
            skip_expccd_chk: removes checks that an image EXPNUM and CCDNUM match those in region file.

        Outputs:
            Makes updates to self.sci.mask in place.
 
        Notes about region files: 
            Example formats for regions:
                fk5;pixel(ra,dec) # expnum=EXPNUM ccdnum=CCDNUM 
                fk5;polgon(ra1,dec1,ra2,dec2,ra3,dec3,....) # expnum=EXPNUM ccdnum=CCDNUM
            Key points:
                "fk5;" is not necessary... coordintes are currrently evaluated assuming they are in the same frame a the image
                "point" and "polygon" are used to distinguish how to further parse information.
                 "(" and ")" denote the portion of the string to be parsed for coordinates
                ra,dec are assumed to be in units of degrees
                comments with format of expnum=VALUE ccdnum=VALUE are assumed to have no spaces and are required unless --skip_expccd_chk is used
        '''

        # kernels for expanding around single-point flags
        xpt_exp=np.array([1,0,-1,1,0,-1,1,0,-1])
        ypt_exp=np.array([1,1,1,0,0,0,-1,-1,-1])

        reg_dict={'nentry':0}
        expnum=self.sci.header['EXPNUM']
        ccdnum=self.sci.header['CCDNUM']
        print("Reading/parsing region mask file: {:s}".format(region_file))
        try:
            rfile=open(region_file,'r')
        except:
            print("Failed to read region file {:s}".format(region_file))
            exit(1)
        # Parse the region file        
        for line in rfile:
            if (line[0] != "#"):
                lbits=line.split()
                expccd_ok=False
                if (skip_expccd_chk):
                    expccd_ok=True
                else:
                    r_exp=None
                    r_ccd=None
                    for bit in lbits:
                        if (re.search("expnum=",bit)):
                            r_exp=int(bit.split("=")[-1])
                        if (re.search("ccdnum=",bit)):
                            r_ccd=int(bit.split("=")[-1])
                    if((r_exp==expnum)and(r_ccd==ccdnum)):
                        expccd_ok=True
                if (expccd_ok):
                    n=reg_dict['nentry']+1
                    reg_dict[n]={}
                    reg_dict['nentry']=n
#                    print("{:s}".format(line))
                    c1=re.search("\(",line).end()
                    c2=re.search("\)",line).start()
                    if (re.search("point",line)):
                        reg_dict[n]['type']="point"
                    if (re.search("polygon",line)):
                        reg_dict[n]['type']="polygon"
                    radec=np.fromstring(line[c1:c2],dtype=float,sep=',')
                    npts=int((radec.size)/2)
                    reg_dict[n]['ra']=np.reshape(radec,(npts,2))[:,0]
                    reg_dict[n]['dec']=np.reshape(radec,(npts,2))[:,1]
                    reg_dict[n]['line']=line
#                else:
#                    print("SKIPPING: {:s}".format(line))
#            else:
#                print("SKIPPING: {:s}".format(line))
        rfile.close()          
      
        # Form set of pixels for each region and add to mask. 
        w = wcsutil.WCS(self.sci.header)
        (ny,nx)=self.sci.mask.shape
#        print(nx,ny)
        print("Identified {:d} regions to be associated with this image".format(reg_dict['nentry']))
        nr_flag_tot=0
        for ireg in range(1,reg_dict['nentry']+1):
#            print(reg_dict[ireg]['line'])
            x, y = w.sky2image(reg_dict[ireg]['ra'],reg_dict[ireg]['dec'])
#            print(x,y)
            if (reg_dict[ireg]['type'] == "point"):
#                print("point")
                xmsk=np.rint(x)+xpt_exp
                ymsk=np.rint(y)+ypt_exp
                xy = np.array([(xmsk[i], ymsk[i]) for i in range(x.size)],dtype=int)
            elif (reg_dict[ireg]['type'] == "polygon"):
#                print("polygon")
                xy = polygon_to_pix(x,y)
#            print(xy)

#           # Check for and remove any portion of the region that would extend beyond the image boundaries
            wsm=np.logical_and(np.logical_and(xy[:,0]-1>=0,xy[:,0]<nx),np.logical_and(xy[:,1]-1>=0,xy[:,1]<ny))
            # Flag remainging pixels    
            self.sci.mask[xy[:,1][wsm]-1,xy[:,0][wsm]-1] |= parse_badpix_mask(reg_flag_val)
            nr_flag=xy[:,1][wsm].size
            nr_flag_tot=nr_flag_tot+nr_flag
            # print("Region flagged {:d} pixels".format(nr_flag))
        print("Finished flagging based on region file.  Total of {:d} (not necessarily unique) pixels flagged as {:s}.".format(nr_flag_tot,reg_flag_val))
        


    @classmethod
    def add_step_args(cls, parser):
        """Add arguments for null_weights and row_zipper
        """
        null_weights.add_step_args(parser)
        row_zipper.add_step_args(parser)
        #parser.add_argument('--custom_weight', action='store_true',default=self.DEFAULT_CUSTOM_WEIGHT,
        #                    help='Run custom weights for STAR and do not write MSK plane for multi-epoch (me)')
        parser.add_argument('--me_wgt_keepmask', action='store', default=cls.DEFAULT_ME_WGT_KEEPMASK,
                            help='Run custom weight for multi-epoch (me) WGT_ME and preserve KEEPMASK')
        parser.add_argument('--null_mask_sci', action='store', default=cls.DEFAULT_NULL_MASK_SCI,
                            help='Names of mask bits to null (or an integer mask) on the SCI plane')
        parser.add_argument('--headfile', action='store', default=cls.DEFAULT_HEADFILE,
                            help='Headfile (containing most update information)')
        parser.add_argument('--hdupcfg', action='store', default=cls.DEFAULT_HDUPCFG,
                            help='Configuration file for header update')
        parser.add_argument('--tilename', action='store', default=cls.DEFAULT_TILENAME,
                            help='Add (optional) TILENAME to SCI header')
        parser.add_argument('--tileid', action='store', type=int, default=cls.DEFAULT_TILEID,
                            help='Add (optional) TILEID to SCI header')

        # Options for the extra streak maskig
        parser.add_argument('--streak_file', action='store', type=str, default='',
                            help='Streak table file path')
        parser.add_argument('--add_width', action='store', type=float, default=0.,
                            help='Broaden streak width by this value (pixels)')
        parser.add_argument('--add_length', action='store', type=float, default=100.,
                            help='Extend streak endpoints by this value (pixels)')
        parser.add_argument('--max_extrapolate', action='store', default=0.0, type=float,
                            help='Do not use streaks extrapolated more than this many degrees')

        # Options to utilize mask regions (from SWarp combine=CLIPPED + clip_mask.py)
        parser.add_argument('--region_mask', action='store', type=str, default='',
                            help='Region file defining masks to be applied to single-epoch images')
        parser.add_argument('--region_skip_imgchk', action='store', type=str, default=cls.DEFAULT_SKIP_IMGCHK,
                            help='Flag to skip checks on expnum/ccdnum when parsing the region_mask')
        parser.add_argument('--reg_flag_val', action='store', type=str, default="CRAY", 
                            help='DESDM BADPIX flag bit to apply to pixels in region_mask (default=CRAY)')


def get_safe_boolean(name, config, config_section):

    """ Get boolean first and if fail, get from config"""
    try:
        param = config.getboolean(config_section, name)
    except:
        param = config.get(config_section, name)
    return param

if __name__ == '__main__':
    CoaddZipperInterpNullWeight.main()
