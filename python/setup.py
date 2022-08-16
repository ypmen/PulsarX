#!/usr/bin/env python
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "filtools",
        sources=["cpp/filtools.cpp",
                 "../src/container/databuffer.cpp",
                 "../src/module/downsample.cpp",
                 "../src/module/preprocesslite.cpp",
                 "../src/module/equalize.cpp",
                 "../src/module/rfi.cpp",
                 "../src/container/AVL.cpp",
                 "../src/utils/utils.cpp"
        ],
        extra_compile_args=["-O3", "-mavx2", "-mfma"],
        library_dirs=["../lib/"],
        libraries=["ymw16", "sofa_c", "boost_log", "fftw3f_threads", "fftw3f", "fftw3_threads", "fftw3"],
        include_dirs=["../include", "../"]
    ),
    Pybind11Extension(
        "orbit_utils",
        sources=["cpp/orbit_utils.cpp"],
	extra_compile_args=["-O3", "-mavx2"],
    ),
]

setup(ext_modules=ext_modules)
