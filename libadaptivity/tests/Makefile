#    Copyright (C) 2009 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Gerard Gorman
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    g.gorman@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA
SHELL = /bin/sh

CPPFLAGS = -I/opt/local/include/vtk-5.10  -DHAVE_VTK=1 -I../load_balance/include \
	-I../adapt3d/include -I../metric_field/include -I../include

F77     = mpif90
FFLAGS  = $(CPPFLAGS) -ffast-math -frecord-marker=4  -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8

FC      = mpif90
FCFLAGS = $(CPPFLAGS) -ffast-math -frecord-marker=4  -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8

CC      = mpicc
CCFLAGS = $(CPPFLAGS) @CCFLAGS@

CXX     = mpicxx
CXXFLAGS= $(CPPFLAGS)   -O3 

LINKER  = mpicxx
LDFLAGS =   
 
LIBS    = ../lib/libadaptivity.a  -llapack -lblas -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkexpat -lvtkzlib -lvtkCommon -ldl -lpthread -lm -lstdc++  -L/opt/local/lib/vtk-5.10  -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4 -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4/../../.. -lgfortran -lquadmath -lm -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4 -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4/../../.. -lgfortran -lquadmath -lm -L./lib -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4 -L/opt/local/lib/gcc47/gcc/x86_64-apple-darwin12/4.7.4/../../.. -lgfortran -lquadmath -lm -lmpi_f77

MAKE    = gmake
INSTALL = /opt/local/bin/ginstall -c

.SUFFIXES: .F .F90 .c .o .a

.F.o:
	$(F77) $(FFLAGS) -c $<

.F90.o:
	$(FC) $(FCFLAGS) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

OBJS = test_adapt_full.o test_metrictensor.o
TESTS = test_adapt_full test_metrictensor

default: $(OBJS)
	$(CXX) $(LDFLAGS) -o test_adapt_full test_adapt_full.o $(LIBS)
	$(CXX) $(LDFLAGS) -o test_metrictensor test_metrictensor.o $(LIBS)

clean:
	rm -f *.o $(TESTS)
