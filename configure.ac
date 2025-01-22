#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([PulsarX], [0.0], [ypmen@mpifr-bonn.mpg.de])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_SRCDIR([include/gridsearch.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
CXXFLAGS='-std=c++11 -O3'

LT_INIT

# Checks for libraries.
AX_BOOST_BASE([1.56], [LIBS="$LIBS -lboost_system -lboost_program_options -lboost_log"], [echo "Boost library was not found!";exit -1])
AC_CHECK_LIB([cfitsio], [ffgky], [], [echo "cfitsio library was not found!";exit -1])
AC_CHECK_LIB([png], [png_set_flush], [], [echo "library not found!";exit -1])
AX_PTHREAD([LIBS+=" -lpthread"], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3], [fftw_malloc], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3f], [fftwf_malloc], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([fftw3f_threads], [fftwf_init_threads], [], [echo "library not found!";exit -1])
AC_CHECK_LIB([lapack], [dgesv_], [], [echo "library not found!";exit -1])

X_LIB_PLOTX([echo "library found"], [echo "library not found!";exit -1])
X_LIB_SOFA([echo "library found"], [echo "library not found!"])

LIBS="$LIBS -lgfortran -lm -lquadmath -lX11"

AX_CHECK_CUDA
AC_OPENMP
AX_CHECK_X86_FEATURES
AX_PYTHON

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MKTIME
AC_CHECK_FUNCS([floor gettimeofday memset strchr])

AC_CONFIG_FILES([Makefile src/Makefile src/pulsar/Makefile])
AC_OUTPUT
