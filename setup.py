#!/usr/bin/env python3 

from numpy.distutils.core import setup, Extension

packages = ['akima86']
package_dir = {'akima86' : 'akima86'}

srcs = ['src/uvipia_omp.f90']
ext = Extension("uvipia_omp", srcs,
                extra_f90_compile_args=['-march=native', '-fopenmp'],
                extra_link_args = ['-lgomp'])
         

setup(name = 'akima86',
        version         = '1.0',
        description     = 'Improved Akima 1-D Interpolation Method',
        author          = 'Hitoshi Akima (original), David Ryglicki (F2003/OpenMP/PY implementation)',
        author_email    = 'david.ryglicki@gmail.com',
        ext_modules     = [ext],
        packages        = packages,
        package_dir     = package_dir,
        install_requires = ['numpy', 'matplotlib'],
        )
