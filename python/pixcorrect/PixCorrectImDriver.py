#!/usr/bin/env python
"""Common code for single steps of pixcorrect-im
"""

# imports
import ctypes
import sys
import logging
from os import path
from ConfigParser import SafeConfigParser 
from argparse import ArgumentParser

import numpy as np

from pixcorrect import corr_util
from pixcorrect.corr_util import logger
from pixcorrect import proddir
from despyfits.DESImage import DESImage, DESImageCStruct

# constants
# exception classes
# interface functions
# classes

class PixCorrectImDriver(object):
    description = None
    step_name = None
    
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
    def parser(cls):
        """Generate a parser for a specific step
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
        parser.add_argument('-i', '--in', nargs=1, 
                                 default=None,
                                 help='input image file name')
        parser.add_argument('-o', '--out', nargs=1, 
                                 default=None,
                                 help='output image file name')

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
                    
        with open(args.saveconfig, 'w') as out_config:
            config.write(out_config)

        return config, args

    @classmethod
    def main(cls):
        config, args = cls.config()

        # start logger
        logging.basicConfig(filename=args.log,
                            format="%(asctime)s\t%(message)s",
                            level=logging.WARNING)
        global logger
        logger = logging.getLogger()
        if args.verbose > 0:
            verbosity = logging.INFO if args.verbose==1 else logging.DEBUG
            logger.setLevel(verbosity)

        logger.addHandler(logging.StreamHandler())
        
        ret_val = cls.run(config)
        
        try:
            ret_val = cls.run(config)
            exit_status = 0 if ret_val is None else ret_val
            sys.exit(ret_val)
        except:
            # If we want to set specific exit status values
            # based on what exceptions get thrown, do that
            # here
            raise

class PixCorrectImStep(PixCorrectImDriver):

    @classmethod
    def __call__(cls, image):
        """Execute the step

        :Parameters:
            - `image`: the DESImage on which to operate

        Applies the correction "in place"
        """
        raise NotImplemetedError

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
        image = DESImage.load(in_fname)

        try:
            if image.mask is None:
                image.init_mask()
        except AttributeError:
            image.init_mask()

        try:
            if image.weight is None:
                image.init_weight()
        except AttributeError:
            image.init_weight()

        ret_code = cls.step_run(image, config)

        out_fname = config.get(cls.step_name, 'out')
        image.save(out_fname)
        
        return ret_code


# internal functions & classes
