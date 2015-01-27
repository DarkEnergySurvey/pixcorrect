"""Run testing code

Runs tests listed in the supplied config file. All classes must be in PYTHONPATH
"""
import unittest
import sys
import os
from doctest import DocTestSuite, DocFileSuite
from ConfigParser import ConfigParser
from argparse import ArgumentParser

sys.path.insert(0, os.path.join(os.path.dirname(__file__),'unittests'))

def main():
    parser = ArgumentParser(description='Run tests.')
    parser.add_argument('config', nargs='?', 
                        default="test.config", help="the configuration file")
    parser.add_argument('--verbose', '-v', action='count', 
                        help="verbosity to be passed to TestRunner")
    args = parser.parse_args()

    config = ConfigParser()
    config.optionxform = str ;# needed to make config args case sensative
    config.read(args.config)

    suite = unittest.TestSuite()

    def true_opts_in_section(section):
        if config.has_section(section):
            return filter(lambda o: config.getboolean(section, o),
                          config.options(section))
        else:
            return []

    # Run unit tests
    for test_name in true_opts_in_section('unittest'):
        test_module = __import__(test_name)
        test_class = getattr(test_module, test_name)
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(test_class))

    # Run doctests embedded in the docstrings in the python modules
    for tested_module in true_opts_in_section('modules'):
        suite.addTests(DocTestSuite(tested_module))

    # Run doctests embedded in reStructuredText text documents
    for tested_doc in true_opts_in_section('docs'):
        doc_full_path = os.path.join(os.path.dirname(__file__), tested_doc)
        suite.addTests(DocFileSuite(doc_full_path))

    unittest.TextTestRunner(verbosity=args.verbose).run(suite)
    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)
