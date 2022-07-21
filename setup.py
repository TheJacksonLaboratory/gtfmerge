#!/usr/bin/env python

import setuptools
import Cython.Build
import sys

"""
run with: python setup.py build_ext --inplace
"""

setuptools.setup(
    ext_modules = [
        Cython.Build.cythonize("gtf.pyx"),
        Cython.Build.cythonize("resolve.pyx"),
    ]
)

