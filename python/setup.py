#!/usr/bin/env python
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "orbit_utils",
        sources=["cpp/orbit_utils.cpp"],
	extra_compile_args=["-O3"],
    ),
]

setup(ext_modules=ext_modules)
