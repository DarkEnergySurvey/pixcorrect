import os
import distutils
from distutils.core import setup, Extension
import glob

bin_files = glob.glob('bin/*')
#inc_files = glob.glob("include/*.h")
#doc_files = glob.glob("doc/*.*") + glob.glob("doc/*/*")


libbiascorrect = Extension(
    'biascorrect',
    sources = ['src/libbiascorrect.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libbpm = Extension(
    'bpm',
    sources = ['src/libbpm.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libfixcol = Extension(
    'fixcol',
    sources = ['src/libfixcol.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libflatcorrect = Extension(
    'flatcorrect',
    sources = ['src/libflatcorrect.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libmasksatr = Extension(
    'masksatr',
    sources = ['src/libmasksatr.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libfpnumber = Extension(
    'fpnumber',
    sources = ['src/libfpnumber.c'],
    include_dirs = ['include', '%s/include' % os.environ['IMSUPPORT_DIR'], '%s/include' % os.environ['DESPYFITS_DIR']],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

# The main call
setup(name='pixcorrect',
      version ='3.1.6',
      description = "Pixel-level image correction",
      author = "Eric Neilsen, Felipe Menanteau",
      author_email = "neilsen@fnal.gov",
      ext_modules = [libbiascorrect, libbpm, libfixcol, libflatcorrect, libmasksatr, libfpnumber],
      packages = ['pixcorrect'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=[ ('ups',['ups/pixcorrect.table']),
                   #('doc', doc_files),
                   #('include', inc_files),
                   ]
      )
