#!/usr/bin/env python

"""
run with 'python setup.py build_ext --inplace'
"""

import setuptools
import Cython.Build

setuptools.setup(
    ext_modules=Cython.Build.cythonize([
        'assign_ids.py',
        'gtf.py',
        'initial2.py',
        'initial1.py',
        'output2.py',
        'output1.py',
        'resolve.py'
    ])
)

