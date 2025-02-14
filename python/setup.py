#!/usr/bin/env python
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "filtools",
        sources=["cpp/filtools.cpp",
        ],
        extra_compile_args=["-O3", "-mavx2", "-mfma"],
        library_dirs=["../../lib/"],
        libraries=["xcontainer", "xmodule", "xutils", "ymw16", "sofa_c", "boost_log", "fftw3f_threads", "fftw3f", "fftw3_threads", "fftw3"],
        include_dirs=["../../include", "../include", "../"]
    ),
    Pybind11Extension(
        "psrdata_reader",
        sources=["cpp/psrdata_reader.cpp",
	],
	extra_compile_args=["-O3", "-mavx2", "-mfma"],
	library_dirs=["../../lib/"],
	libraries=["xcontainer", "xformats", "xmodule", "xutils", "cfitsio", "ymw16", "sofa_c", "boost_log", "fftw3f_threads", "fftw3f", "fftw3_threads", "fftw3"],
	include_dirs=["../../include", "../include", "../"]
    ),
    Pybind11Extension(
        "orbit_utils",
        sources=["cpp/orbit_utils.cpp", "cpp/numerical_algorithm.cpp"],
		extra_compile_args=["-O3", "-mavx2"],
		libraries=["lapack"],
        include_dirs=["../../include", "../include", "../"]
    ),
    Pybind11Extension(
        "formats",
        sources=["cpp/formats.cpp",
	],
        extra_compile_args=["-O3", "-mavx2", "-mfma"],
	library_dirs=["../../lib/"],
	libraries=["xformats", "cfitsio"],
	include_dirs=["../../include", "../include", "../"]
    ),
]

setup(ext_modules=ext_modules)
