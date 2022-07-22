#!/usr/bin/env python

"""
run with 'python setup.py build_ext --inplace'
"""

import Cython.Build
import setuptools

setuptools.setup(
    ext_modules=Cython.Build.cythonize([
        'assign_ids.pyx',
        'gtf.pyx',
        'initial2.pyx',
        'initial.pyx',
        'output2.pyx',
        'output.pyx',
        'resolve.pyx',  
    ])
)

