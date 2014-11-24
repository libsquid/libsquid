#!/usr/bin/env python

from distutils.core import setup, Extension

my_srcs=['pysquid.i',
         '../libsquid_base.c',
         '../libsquid_projections.c',
         '../libsquid_utils.c']
pysquid_module = Extension('_pysquid',
                           sources=my_srcs,
                           include_dirs=['../'])

setup(name='pysquid',
      author='Jim Wren',
      description='LibSQUID Library Interface for Python',
      version='0.5',
      ext_modules=[pysquid_module],
      py_modules=["pysquid"],
)
                                        
