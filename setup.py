#!/usr/bin/env python

import setuptools
import Cython.Build 
import sys

"""
run with: python setup.py build_ext --inplace
"""

setuptools.setup(
    ext_modules=Cython.Build.cythonize([
        "assign_ids.pyx",
        "gtf.pyx", 
        "initial.pyx",
        "initial2.pyx",
        "output.pyx", 
        "output2.pyx",
        "resolve.pyx"]
    )
)

