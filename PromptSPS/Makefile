# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
# You can find plugins here: INSERT URL
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files
THEPEGINCLUDE=-I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/thepeg/1.9.2p1-ikhhed3/include
GSLINCLUDE=-I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/1.16/include/
HERWIGINCLUDE=-I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/herwigpp/2.7.1-ikhhed3/include

LDFLAGS= 
SHARED_FLAG=-shared 
INCLUDE = $(THEPEGINCLUDE) $(GSLINCLUDE) $(HERWIGINCLUDE)
#
# C++ flags
# 
CXX=g++
CXXFLAGS=-O2 -fuse-cxa-atexit -g

ALLCCFILES=$(shell echo *.cc)

default : MEgg2JpsiJpsi.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c -shared $< -o $@  

MEgg2JpsiJpsi.so : MEgg2JpsiJpsi.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) $(SHARED_FLAG) $(LDFLAGS) MEgg2JpsiJpsi.o \
	 -o MEgg2JpsiJpsi.so


clean:
	rm -f $(ALLCCFILES:.cc=.so)  $(ALLCCFILES:.cc=.o) 
