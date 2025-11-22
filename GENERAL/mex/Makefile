#
# SDPA-M: $Revision: 7.3 $
#
# Linux: using MATLAB built-in BLAS and LAPACK

-include ../etc/make.inc
OS = $(shell uname)

# Set manually the following MEX full path
MEX = "auto"
# for Linux (depends on version)
#MEX = /usr/local/matlab/MatlabR2008b/bin/mex
# for Mac (depends on version)
ifeq ($(OS),Darwin)
	MEX = /Applications/MATLAB_R2008b.app/bin/mex
endif


# Automatic search MEX full path by invoking matlab
ifeq ($(MEX),"auto")
	MEX = ${shell matlab -nodisplay -r "s=matlabroot; fprintf('%s\n',s);exit" | tail -2 | head -1 }/bin/mex
endif

MATLAB_LAPACK = -lmwlapack -lmwblas
DEF_FPIC = 
ARCH = $(shell uname -m)
ifeq ($(ARCH),x86_64)
        DEF_FPIC = -fPIC
endif

OPTIONS = CXXFLAGS="-Wall  ${DEF_FPIC}" -O -DNDEBUG CC=${CC} CXX=${CXX}
ALL_INCLUDE = -I${SDPA_DIR}/include  -I${MUMPS_DIR}/include 
ALL_LIBS    = ${SDPA_LIB} ${MUMPS_LIBS} ${MATLAB_LAPACK} ${FLIBS}

all: mexBinaryCheck mexsdpa mexSedumiWrap mexAggSDPcones mexDisAggSDPsol
# all: mexSedumiWrap

mexBinaryCheck:
	@echo "---------------------------------------------"
	@echo "Mex Compilation Start                        "
	@echo "---------------------------------------------"
	@if [ -f $(MEX) ]; then \
	echo "---------------------------------------------" ;\
	echo " We use                                      " ;\
	echo "   " $(MEX) ;\
	echo "   for Mex compilaton                        " ;\
	echo "---------------------------------------------" ;\
	else \
	echo "---------------------------------------------" ;\
	echo "Automatic search is failed,                  " ;\
	echo " Set the full path for mex command mannually " ;\
	echo " in "`pwd`"/Makekefile" ;\
	echo " and try 'make' in `pwd`                     " ;\
	echo "---------------------------------------------" ;\
	false ;\
	fi

mexsdpa:
	${MEX} -largeArrayDims \
	  ${OPTIONS} ${ALL_INCLUDE} \
	  mexsdpa.cpp mexFprintf.c \
	  ${ALL_LIBS}

mexSedumiWrap:
	${MEX} -largeArrayDims \
	  ${OPTIONS} ${ALL_INCLUDE} \
	  mexSedumiWrap.cpp mexFprintf.c \
	  ${ALL_LIBS}

mexAggSDPcones:
	${MEX} -largeArrayDims \
	  ${OPTIONS} \
	  mexAggSDPcones.cpp mexFprintf.c

mexDisAggSDPsol:
	${MEX} -largeArrayDims \
	  ${OPTIONS} \
	  mexDisAggSDPsol.cpp mexFprintf.c


# clean all generated files
clean:
	rm -f *~
	rm -f mexsdpa.mex* mexSedumiWrap.mex*
	rm -f mexAggSDPcones.mex* mexDisAggSDPsol.mex*
	rm -f mexFprintf.mex*

distclean: clean

#
# End of File
#
