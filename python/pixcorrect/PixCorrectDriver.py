#!/usr/bin/env python
"""Common code for single steps of pixcorrect-im
"""

# imports
import ctypes
import sys
import logging

from os import path
from ConfigParser import SafeConfigParser, NoOptionError
from argparse import ArgumentParser

import numpy as np

from pixcorrect import corr_util
from pixcorrect.corr_util import logger
from pixcorrect import proddir
from despyfits.DESImage import DESImage, DESBPMImage, DESImageCStruct
from despyfits.DESFocalPlaneImages import DESFocalPlaneImages

# constants
# exception classes
# interface functions
# classes

class PixCorrectDriver(object):

    @classmethod
    def run(cls, config):
        """Customized execution for this step

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """
        raise NotImplemetedError

    @classmethod
    def add_step_args(cls, parser):
        """Add arguments specific to this step
        """
        pass


    @classmethod
    def common_parser(cls):
        """Generate a parser for a specific step
        """
        default_config = path.join(proddir, 'etc', cls.step_name+'.config')
        # default_out_config = path.join(cls.step_name+'-as_run'+'.config')
        default_out_config = ""

        # Argument parser
        parser = ArgumentParser(description=cls.description)
        parser.add_argument("config", default=default_config, nargs="?",
                            help="Configuration file filename")
        parser.add_argument('-s', '--saveconfig', 
                                 default=default_out_config,
                                 help="output config file")
        parser.add_argument('-l', '--log', 
                                 default="", 
                                 help="the name of the logfile")
        parser.add_argument('-v', '--verbose', action="count", 
                                 help="be verbose")

        return parser

    @classmethod
    def parser(cls):
        """Generate a parser for a specific step
        """
        parser = cls.common_parser()
        cls.add_step_args(parser)

        return parser


    @classmethod
    def config(cls):
        """Return a configuration object for the step
        """
        args = cls.parser().parse_args()
        
        # load configuration
        config = SafeConfigParser() 
        config.read(args.config) 

        section = cls.step_name
        if not config.has_section(section):
            config.add_section(section)

        for argument, value in args._get_kwargs():
            value = getattr(args, argument)
            if value is not None:
                if type(value)==type([]):
                    value=value[0]
                config.set(section, argument, str(value))
                    
        if args.saveconfig is not None and len(args.saveconfig)>0:
            with open(args.saveconfig, 'w') as out_config:
                config.write(out_config)

        return config, args

    @classmethod
    def main(cls):
        config, args = cls.config()

        # start logger
        if args.log is not None and len(args.log)>0:
            logging.basicConfig(filename=args.log,
                                format="%(asctime)s %(levelname)s:\t%(message)s",
                                level=logging.WARNING)
            sh = logging.StreamHandler()
            sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s:\t%(message)s"))
            logger.addHandler(sh)
        else:
            logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                                level=logging.WARNING)

        global logger
        logger = logging.getLogger()
        if args.verbose > 0:
            verbosity = logging.INFO if args.verbose==1 else logging.DEBUG
            logger.setLevel(verbosity)

        
        try:
            ret_val = cls.run(config)
            exit_status = 0 if ret_val is None else ret_val
            sys.exit(ret_val)
        except:
            # If we want to set specific exit status values
            # based on what exceptions get thrown, do that
            # here
            raise

class PixCorrectStep(PixCorrectDriver):

    # In all subclasses, the __call__ method should take all unpacked
    # parameters, while the step_run takes the sci image (or array of sci images)
    # and the configuration object

    # The call method should actually apply the corrections
    @classmethod
    def __call__(cls):
        """Execute the step

        :Parameters:
            - `image`: the DESImage on which to operate

        Applies the correction "in place"
        """
        raise NotImplemetedError


    # The step_run method unpacks parameters from config, and 
    # calls __call__ to do the corrections.
    @classmethod
    def step_run(cls, config):
        """Customized execution for this step

        :Parameters:
            - `config`: the configuration from which to get other parameters

        """
        ret_code = cls.__call__()
        return ret_code

    @classmethod
    def run(cls, config):
        """Execute the step, loading and running the input and output
        """
        ret_code = cls.step_run(config)
        return ret_code

class PixCorrectImStep(PixCorrectStep):

    @classmethod
    def parser(cls):
        """Generate a parser for a specific step
        """
        parser = cls.common_parser()

        parser.add_argument('-i', '--in', 
                                 default=None,
                                 help='input image file name')
        parser.add_argument('-o', '--out', 
                                 default=None,
                                 help='output image file name')
        """
                            parser.add_argument('-n', '--ccdnum', nargs='?', 
                            type=int,
                            help='input image CCD number')
        """
        cls.add_step_args(parser)

        return parser

    # The step_run method unpacks parameters from config, and 
    # calls __call__ to do the corrections.
    @classmethod
    def step_run(cls, image, config):
        """Customized execution for this step

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        ret_code = cls.__call__(image)
        return ret_code

    @classmethod
    def run(cls, config):
        """Execute the step, loading and running the input and output
        """
        in_fname = config.get(cls.step_name, 'in')
        try:
            ccdnum = config.getint(cls.step_name, 'ccdnum')
            image = DESImage.load(in_fname, ccdnum=ccdnum)
        except NoOptionError:
            image = DESImage.load(in_fname)

            
        ret_code = cls.step_run(image, config)

        out_fname = config.get(cls.step_name, 'out')
        image.save(out_fname)
        
        return ret_code

class PixCorrectFPStep(PixCorrectStep):

    @classmethod
    def parser(cls):
        """Generate a parser for a specific step
        """
        parser = cls.common_parser()

        parser.add_argument('-i', '--in', 
                                 default=None,
                                 help='input image file name')
        parser.add_argument('-o', '--out', 
                                 default=None,
                                 help='output image file name')

        cls.add_step_args(parser)

        return parser

    # The step_run method unpacks parameters from config, and 
    # calls __call__ to do the corrections.
    @classmethod
    def step_run(cls, image, config):
        """Customized execution for this step

        :Parameters:
            - `image`: the DESImage on which to operate
            - `config`: the configuration from which to get other parameters

        """
        ret_code = cls.__call__(image)
        return ret_code


    @classmethod
    def run(cls, config):
        """Execute the step, loading and running the input and output
        """
        in_fname = config.get(cls.step_name, 'in')
        images = DESFocalPlaneImages.load(in_fname)

        ret_code = cls.step_run(images, config)

        out_fname_template = config.get(cls.step_name, 'out')
        images.save(out_fname_template)
        
        return ret_code

class PixCorrectMultistep(PixCorrectDriver):

    def __init__(self, config):
        self.config = config
        self._image_data = {}

    @classmethod
    def parser(cls):
        """Generate a parser for a specific step
        """
        parser = cls.common_parser()

        parser.add_argument('-i', '--in', 
                                 default=None,
                                 help='input image file name')
        parser.add_argument('-o', '--out', 
                                 default=None,
                                 help='output image file name')

        cls.add_step_args(parser)

        return parser

    @classmethod
    def run(cls, config):
        config.set(cls.config_section, 'sci', 
                   config.get(cls.config_section, 'in'))
        pix_corrector = cls(config)
        ret_value = pix_corrector()
        return ret_value

    def image_data(self, fname):
        raise NotImplemetedError

    def __getattr__(self, image_name):
        """Create a shortcut to images using object attributes
        """
        return self.image_data(image_name)

    def clean_im(self, image_name):
        """Let python garbage collect the memory used for an image

        :Parameters:
            -`image_name`: the type of image to clean
        """
        if image_name in self._image_data:
            del self._image_data[image_name]
        

    def do_step(self, step_name):
        if not self.config.has_option(self.config_section, step_name):
            return False

        try:
            # If the parameter is a boolean, interpret is
            # as an on/off switch
            doit = self.config.getboolean(self.config_section, step_name)
            return doit
        except:
            # Otherwise, interpret it as a value associated with
            # the step, and assume we want to perform the step
            return True



# internal functions & classes

#
# RAG: new function to read a file (containing a list of files) and tender the results as a python list.
# RAG: likely this should be moved to despymisc (have left it here for testing).
#
def filelist_to_list(input_file_list, column_used=0, delimeter=None, check_files_exist=True, append_missing_files=False):
    """Create a list of file names from a file containing a list of file names.
        :Parameters:
            -`input_file_list`: the filename containing the list of files
            -`column_used': column within file from which to draw name(s) (default=0)
            -`delimeter': delimeter for parsing columns (default=None --> whitespace)
            -`check_files_exist': function will check that files exist before adding to list (default=True)
            -`append_missing_files': function appends files even if they are missing (default=False)
        :Returns: list 
    """

    list_of_files=[]
    try:
        f_listfile=open(input_file_list,'r')
    except:
        raise IOError("File not found.  Missing input list %s " % input_file_list)
    for line in f_listfile:
        line=line.strip()
        if (delimeter is None):
            columns=line.split()
        else:
            columns=line.split(delimeter)
        FileExists=True
        if (check_files_exist):
            if (not(path.isfile(columns[column_used]))):
                FileExists=False
        if ((append_missing_files)or(FileExists)):
            list_of_files.append(columns[column_used])
    f_listfile.close()
    return(list_of_files)

