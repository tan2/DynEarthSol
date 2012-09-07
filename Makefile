# -*- Makefile -*-
#
# Makefile for DynEarthSol3D
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#

## Execute "make" if making production run. Or "make debug=1" for debugging run.
##
## dim = 3: 3D code; 2: 2D code
## debug = 0: optimized build; 1: debugging build

dim = 3
debug = 0

## Select C++ compiler
CXX = g++


########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

CXXFLAGS = -g
LDFLAGS = -lm

ifeq ($(debug), 0)
	CXXFLAGS += -O2
else
	CXXFLAGS += -O0 -Wall
endif

##

SRCS =	\
	dynearthsol3d.cxx


INCS =	\
	constants.h \
	parameters.h \


OBJS = $(SRCS:.cxx=.o)

EXE = dynearthsol3d


## Libraries

TET_SRCS = tetgen/predicates.cxx tetgen/tetgen.cxx
TET_INCS = tetgen/tetgen.h
TET_OBJS = $(TET_SRCS:.cxx=.o)

TRI_SRCS = triangle/triangle.c
TRI_INCS = triangle/triangle.h
TRI_OBJS = $(TRI_SRCS:.c=.o)

## Action

all: $(EXE)

$(EXE): $(TET_OBJS) $(TRI_OBJS) $(OBJS)
	$(CXX) $(LDFLAGS) $(TRI_OBJS) $(OBJS) -o $@
	@# snapshot of the code for building the executable
	@which hg 2>&1 > /dev/null && (hg summary; hg diff) > snapshot.diff

$(OBJS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $<

$(TRI_OBJS): %.o : %.c $(TRI_INCS)
	@# Triangle cannot be compiled with -O2
	$(CXX) $(CXXFLAGS) -O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS -c $< -o $@

tetgen/predicates.o: tetgen/predicates.cxx $(TET_INCS)
	@# Compiling J. Shewchuk predicates, should always be
	@# equal to -O0 (no optimization). Otherwise, TetGen may not
	@# work properly.
	$(CXX) $(CXXFLAGS) -DTETLIBRARY -O0 -c $< -o $@

tetgen/tetgen.o: tetgen/tetgen.cxx $(TET_INCS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DTETLIBRARY -O1 -c $< -o $@

clean:
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
