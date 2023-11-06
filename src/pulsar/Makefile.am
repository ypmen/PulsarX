bin_PROGRAMS=dedisperse_all_fil psrfold_fil psrfold_fil2 filtool dmffdot chop8 pulsar_injector

AM_CPPFLAGS=-I$(top_srcdir)/include
LDADD=-lboost_program_options -lxmodule -lxutils -lxcontainer -lxformats -lymw16

dedisperse_all_fil_SOURCES=dedisperse_all_fil.cpp pulsarsearch.cpp
psrfold_fil_SOURCES=predictor.cpp dedispersionliteU.cpp archivelite.cpp gridsearch.cpp archivewriter.cpp psrfold_fil.cpp
psrfold_fil2_SOURCES=predictor.cpp archivelite.cpp gridsearch.cpp archivewriter.cpp psrfold_fil2.cpp
filtool_SOURCES=filtool.cpp filmaker.cpp
dmffdot_SOURCES=predictor.cpp dmffdot.cpp archivelite.cpp gridsearch.cpp pulsarplot.cpp archivewriter.cpp
chop8_SOURCES=chop8.cpp
pulsar_injector_SOURCES=pulsar_injector.cpp predictor.cpp

psrfold_fil_SOURCES+=pulsarplot.cpp
psrfold_fil2_SOURCES+=pulsarplot.cpp
