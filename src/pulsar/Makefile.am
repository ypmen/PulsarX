bin_PROGRAMS=dedisperse_all_fil psrfold_fil psrfold_fil2 filtool dmffdot

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module -L$(top_srcdir)/src/ymw16
LDADD=-lboost_program_options -lmodule -lutils -lcontainer -lformats -lymw16

dedisperse_all_fil_SOURCES=dedisperse_all_fil.cpp pulsarsearch.cpp
psrfold_fil_SOURCES=predictor.cpp dedispersionliteU.cpp archivelite.cpp gridsearch.cpp psrfold_fil.cpp
psrfold_fil2_SOURCES=predictor.cpp dedispersionX.cpp archivelite.cpp gridsearch.cpp psrfold_fil2.cpp
filtool_SOURCES=filtool.cpp filmaker.cpp
dmffdot_SOURCES=predictor.cpp dmffdot.cpp archivelite.cpp gridsearch.cpp pulsarplot.cpp

psrfold_fil_SOURCES+=pulsarplot.cpp
psrfold_fil2_SOURCES+=pulsarplot.cpp
