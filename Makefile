# -*- Makefile -*-
#
# Makefile for DynEarthSol3D
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#

## Execute "make" if making production run. Or "make opt=0 openmp=0" for debugging run.
##
## ndims = 3: 3D code; 2: 2D code
## opt = 1 ~ 3: optimized build; others: debugging build
## openmp = 1: enable OpenMP
## gprof = 1: enable instrumentation for GNU profiler

ndims = 2
opt = 2
openmp = 1
gprof = 0

## Select C++ compiler
CXX = g++

## Boost location and library name
BOOST_ROOT_DIR =
BOOST_LDFLAGS = -lboost_program_options

########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

ifdef BOOST_ROOT_DIR
	BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)/include
	BOOST_LDFLAGS += -L$(BOOST_ROOT_DIR)/lib -Wl,--rpath=$(BOOST_ROOT_DIR)/lib
endif

ifeq ($(CXX), g++)
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS = -std=c++0x -march=native ## -g is incompatible with -flto
		CXXFLAGS += -Ofast -funroll-loops -fwhole-program -flto
		LDFLAGS += -flto
	else
		CXXFLAGS += -O0 -Wall -Wno-unused-function -Wno-unknown-pragmas -fbounds-check -ftrapv
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp -DUSE_OMP
		LDFLAGS += -fopenmp
	endif

	ifeq ($(gprof), 1)
		CXXFLAGS += -pg
		LDFLAGS += -pg
	endif
else
# the only way to display the error message in Makefile ...
all:
	@echo "Unknown compiler, check the definition of 'CXX' in the Makefile."
	@false
endif

##

SRCS =	\
	barycentric-fn.cxx \
	brc-interpolation.cxx \
	bc.cxx \
	dynearthsol.cxx \
	fields.cxx \
	geometry.cxx \
	ic.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	nn-interpolation.cxx \
	output.cxx \
	remeshing.cxx \
	rheology.cxx \
	markerset.cxx

INCS =	\
	array2d.hpp \
	barycentric-fn.hpp \
	constants.hpp \
	parameters.hpp \
	matprops.hpp \
	sortindex.hpp \
	utils.hpp \
	mesh.hpp \
	markerset.hpp

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

C3X3_DIR = 3x3-C
C3X3_LIBNAME = 3x3

ANN_DIR = ann
ANN_LIBNAME = ANN
CXXFLAGS += -I$(ANN_DIR)/include

## Action

all: $(EXE)

$(EXE): $(M_OBJS) $(OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a $(ANN_DIR)/lib/lib$(ANN_LIBNAME).a
	$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) \
		-L$(C3X3_DIR) -l$(C3X3_LIBNAME) -L$(ANN_DIR)/lib -l$(ANN_LIBNAME) -o $@
	@# snapshot of the code for building the executable
	@which hg 2>&1 > /dev/null && (hg summary; hg diff) > snapshot.diff

$(OBJS): %.$(ndims)d.o : %.cxx $(INCS)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -c $< -o $@

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

$(C3X3_DIR)/lib$(C3X3_LIBNAME).a:
	@+$(MAKE) -C $(C3X3_DIR)

$(ANN_DIR)/lib/lib$(ANN_LIBNAME).a:
	@+$(MAKE) -C $(ANN_DIR) linux-g++

deepclean:
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
	@+$(MAKE) -C $(C3X3_DIR) clean

clean:
	@rm -f $(OBJS) $(EXE)

