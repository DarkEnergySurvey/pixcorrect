Automated testing of pixcorrect
===============================

Tests can be written in either of python's two native testing
mechanisms: unittest or doctest.

Tests of both types are run by $PIXCORRECT_DIR/test/test.py. A
configuration file, an example of which can be found in
$PIXCORRECT_DIR/test/test.config, defines which tests are to be run.

It can be run thus::

  $ cd $PIXCORRECT_DIR/test
  $ python test.py test.config


To add unittest tests, the unittest modules should be placed in
$PIXCORRECT_DIR/test, and corresponding entries added to the test.config
files.

If modules have internal documentation that can be tested automaticall
by doctest, the modules can be added to $PIXCORRECT_DIR/test.config
directly.

If there are test documents that can be run using doctest, they should
be added to either $PIXCORRECT_DIR/doctests (if they are stand-alone
testing documents), or $PIXCORRECT_DIR/doc/examples (if they are to be
included as examples in documentation).
