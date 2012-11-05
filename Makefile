# -*- Makefile -*-
#
# Makefile for DynEarthSol3D
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#

## Execute "make" if making production run. Or "make debug=1" for debugging run.
##
## ndims = 3: 3D code; 2: 2D code
## debug = 0: optimized build; 1: debugging build

ndims = 2
debug = 0

## Select C++ compiler
CXX = g++

## Boost
BOOSTCXXFLAGS =
BOOSTLDFLAGS = -lboost_program_options-mt

########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

CXXFLAGS = -g -std=c++0x
LDFLAGS = -lm

ifeq ($(debug), 0)
	CXXFLAGS += -O2 -DBOOST_DISABLE_ASSERTS -DNDEBUG
else
	CXXFLAGS += -O0 -Wall -Wno-unused-function
endif

##

SRCS =	\
	dynearthsol.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	sortindex.cxx

INCS =	\
	constants.hpp \
	parameters.hpp \
	utils.hpp

OBJS = $(SRCS:.cxx=.$(ndims)d.o)

EXE = dynearthsol$(ndims)d


## Libraries

TET_SRCS = tetgen/predicates.cxx tetgen/tetgen.cxx
TET_INCS = tetgen/tetgen.h
TET_OBJS = $(TET_SRCS:.cxx=.o)

TRI_SRCS = triangle/triangle.c
TRI_INCS = triangle/triangle.h
TRI_OBJS = $(TRI_SRCS:.c=.o)

ifeq ($(ndims), 2)
	M_SRCS = $(TRI_SRCS)
	M_INCS = $(TRI_INCS)
	M_OBJS = $(TRI_OBJS)
else
	M_SRCS = $(TET_SRCS)
	M_INCS = $(TET_INCS)
	M_OBJS = $(TET_OBJS)
	CXXFLAGS += -DTHREED
endif

## Action

all: $(EXE)

$(EXE): $(M_OBJS) $(OBJS)
	$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOSTLDFLAGS) -o $@
	@# snapshot of the code for building the executable
	@which hg 2>&1 > /dev/null && (hg summary; hg diff) > snapshot.diff

$(OBJS): %.$(ndims)d.o : %.cxx $(INCS)
	$(CXX) $(CXXFLAGS) $(BOOSTCXXFLAGS) -c $< -o $@

$(TRI_OBJS): %.o : %.c $(TRI_INCS)
	@# Triangle cannot be compiled with -O2
	$(CXX) $(CXXFLAGS) -O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS -c $< -o $@

tetgen/predicates.o: tetgen/predicates.cxx $(TET_INCS)
	@# Compiling J. Shewchuk predicates, should always be
	@# equal to -O0 (no optimization). Otherwise, TetGen may not
	@# work properly.
	$(CXX) $(CXXFLAGS) -DTETLIBRARY -O0 -c $< -o $@

tetgen/tetgen.o: tetgen/tetgen.cxx $(TET_INCS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DTETLIBRARY -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -c $< -o $@

clean:
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
