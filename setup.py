import distutils
from distutils.core import setup
import glob

import shlib 
from shlib.build_shlib import SharedLibrary

bin_files = glob.glob('bin/*')

libbpm = SharedLibrary(
    'bpm',
    sources = ['src/libbpm.c'],
    include_dirs = ['include'],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libfixcol = SharedLibrary(
    'fixcol',
    sources = ['src/libfixcol.c'],
    include_dirs = ['include'],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

libmasksatr = SharedLibrary(
    'masksatr',
    sources = ['src/libmasksatr.c'],
    include_dirs = ['include'],
    extra_compile_args = ['-O3','-g','-Wall','-shared','-fPIC'])

# The main call
setup(name='pixcorrect',
      version ='0.1.1',
      license = "GPL",
      description = "Pixel-level image correction",
      author = "Eric Neilsen",
      author_email = "neilsen@fnal.gov",
      shlibs = [libbpm, libfixcol, libmasksatr],
      packages = ['pixcorrect'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=[('ups',['ups/pixcorrect.table'])],
      )

