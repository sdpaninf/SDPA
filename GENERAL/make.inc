#default environments

CC       = gcc
CXX      = g++
CFLAGS   = -O3
CXXFLAGS = $(CFLAGS)

LAPACK_LIBS = -llapack
BLAS_LIBS   = -lblas

MUMPS_DIR = /usr/local/mumps

# APPEND_LIBS = -lgfortran
# ADDITIONAL_CLEAN = 

# F77_FUNC 
DEF_F77_FUNC = -DF77_FUNC\(name,NAME\)=name\ \#\#\ _ 
# DEF_F77_FUNC = -DF77_FUNC\(name,NAME\)=name\ \#\# 
# DEF_F77_FUNC = -DF77_FUNC\(name,NAME\)=name\ \#\#\ __ 

DEF_FPIC = 

USERNAME = $(shell whoami)
HOSTNAME = $(shell hostname | cut -f1 -d'.')
ARCH     = $(shell uname -m)
ifeq ($(ARCH),x86_64)
	DEF_FPIC = -fPIC
endif

# -----------------------------------------
# environments for each user
# -----------------------------------------

ifeq ($(USERNAME),makoto)
	CFLAGS      = -O3 -pg -g -Wall
	LAPACK_DIR  = $(HOME)/lib/b/lapack
	LAPACK_LIBS = -L$(LAPACK_DIR)/lib -llapack
	BLAS_LIBS   = -L$(LAPACK_DIR)/lib -lf77blas -lcblas -latlas
	MUMPS_DIR   = $(HOME)/lib/b/mumps/single
	APPEND_LIBS = -lgfortran -lpthread
	ADDITIONAL_CLEAN = r r2 gmon.out tmp.dat-s tmp.ini-s
endif
