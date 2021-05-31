bin_PROGRAMS=dedisperse_all dedisperse_all_fil psrfold psrfold_fil filtool

AM_CPPFLAGS=-I$(top_srcdir)/include
LDFLAGS=-L$(top_srcdir)/src/container -L$(top_srcdir)/src/formats -L$(top_srcdir)/src/utils -L$(top_srcdir)/src/module -L$(top_srcdir)/src/ymw16
LDADD=-lboost_program_options -lmodule -lutils -lcontainer -lformats -lymw16

dedisperse_all_SOURCES=dedisperse_all.cpp pulsarsearch.cpp
dedisperse_all_fil_SOURCES=dedisperse_all_fil.cpp pulsarsearch.cpp
psrfold_SOURCES=dedispersionliteU.cpp archivelite.cpp gridsearch.cpp psrfold.cpp
psrfold_fil_SOURCES=dedispersionliteU.cpp archivelite.cpp gridsearch.cpp psrfold_fil.cpp
filtool_SOURCES=filtool.cpp filmaker.cpp

psrfold_SOURCES+=pulsarplot.cpp
psrfold_fil_SOURCES+=pulsarplot.cpp
