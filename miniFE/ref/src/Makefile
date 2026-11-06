#-----------------------------------------------------------------------

MINIFE_TYPES =  \
        -DMINIFE_SCALAR=double   \
        -DMINIFE_LOCAL_ORDINAL=int      \
        -DMINIFE_GLOBAL_ORDINAL=int

#MINIFE_MATRIX_TYPE = -DMINIFE_CSR_MATRIX
MINIFE_MATRIX_TYPE = -DMINIFE_ELL_MATRIX

#-----------------------------------------------------------------------

CFLAGS = -O3 -g -fopt-info-vec-optimized
CXXFLAGS = $(CFLAGS)

CPPFLAGS = -I. -I../utils -I../fem $(MINIFE_TYPES) $(MINIFE_MATRIX_TYPE) \
	-DHAVE_MPI -DMPICH_IGNORE_CXX_SEEK \
	-DMINIFE_REPORT_RUSAGE

LDFLAGS=
LIBS=

CXX=mpicxx
CCi=mpicc

#CXX=g++
#CC=gcc

include make_targets

