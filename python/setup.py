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
        "psrdata_reader",
        sources=["cpp/psrdata_reader.cpp",
		 "../src/formats/filterbank.cpp",
		 "../src/formats/integration.cpp",
		 "../src/formats/hdu.cpp",
		 "../src/formats/psrfits.cpp",
		 "../src/module/filterbankreader.cpp",
		 "../src/container/databuffer.cpp",
		 "../src/module/filterbankwriter.cpp",
		 "../src/module/psrfitsreader.cpp",
		 "../src/container/AVL.cpp",
		 "../src/utils/utils.cpp"
	],
	extra_compile_args=["-O3", "-mavx2", "-mfma"],
	library_dirs=["../lib/"],
	libraries=["cfitsio", "ymw16", "sofa_c", "boost_log", "fftw3f_threads", "fftw3f", "fftw3_threads", "fftw3"],
	include_dirs=["../include", "../"]
    ),
    Pybind11Extension(
        "orbit_utils",
        sources=["cpp/orbit_utils.cpp"],
	extra_compile_args=["-O3", "-mavx2"],
        include_dirs=["../include", "../"]
    ),
    Pybind11Extension(
        "formats",
        sources=["cpp/formats.cpp",
		 "../src/formats/integration.cpp",
		 "../src/formats/hdu.cpp",
		 "../src/formats/psrfits.cpp",
		 "../src/formats/filterbank.cpp"
	],
        extra_compile_args=["-O3", "-mavx2", "-mfma"],
	library_dirs=[],
	libraries=["cfitsio"],
	include_dirs=["../include", "../"]
    ),
]

setup(ext_modules=ext_modules)
