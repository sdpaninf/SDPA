
SDPA_MAKEFILE = ../Makefile
# This file is automatically made by './configure' script
###################################################################
# The following variables are set from 
# the arguments given to './configure' script for building SDPA
###################################################################
SDPA_PACKAGE= ${shell grep "^PACKAGE_STRING =" ${SDPA_MAKEFILE} | awk '{printf("%s.%s",$$3,$$4)}'}

SDPA_DIR    = ${shell grep "^prefix ="       ${SDPA_MAKEFILE} | cut -f2 -d= | sed 's+ ++'}/${SDPA_PACKAGE}
SDPA_LIB    = ${SDPA_DIR}/lib/libsdpa.a 
MUMPS_DIR   = ${SDPA_DIR}/mumps/build
MUMPS_LIBS  = -L${MUMPS_DIR}/lib -ldmumps -lmumps_common -lpord \
              -L${MUMPS_DIR}/libseq -lmpiseq
CC          = ${shell grep "^CC ="          ${SDPA_MAKEFILE} | cut -f2 -d=}
CFLAGS      = ${shell grep "^CFLAGS = "    ${SDPA_MAKEFILE} | cut -f2 -d=}
CXX         = ${shell grep "^CXX ="          ${SDPA_MAKEFILE} | cut -f2 -d=}
CXXFLAGS    = ${shell grep "^CXXFLAGS = "    ${SDPA_MAKEFILE} | cut -f2 -d=}
LAPACK_LIBS = ${shell grep "^LAPACK_LIBS = " ${SDPA_MAKEFILE} | cut -f2 -d=}
BLAS_LIBS   = ${shell grep "^BLAS_LIBS = "   ${SDPA_MAKEFILE} | cut -f2 -d=}
pthread_lib = ${shell grep "^pthread_lib = " ${SDPA_MAKEFILE} | cut -f2 -d=}
FLIBS       = ${shell grep "^FLIBS = "       ${SDPA_MAKEFILE} | cut -f2 -d=}
###################################################################

TARGET = make.inc
all: ${TARGET}
${TARGET}:
	rm -f make.inc
	echo "SDPA_DIR    = " ${SDPA_DIR}  >> ${TARGET}
	echo "SDPA_LIB    = " ${SDPA_LIB}  >> ${TARGET} 
	echo "MUMPS_DIR   = " ${MUMPS_DIR} >> ${TARGET}
	echo "CC          = " ${CC} >> ${TARGET}
	echo "CFLAGS      = " ${CFLAGS} >> ${TARGET}
	echo "CXX         = " ${CXX} >> ${TARGET}
	echo "CXXFLAGS    = " ${CXXFLAGS} >> ${TARGET}
	echo "LAPACK_LIBS = " ${LAPACK_LIBS} >> ${TARGET}
	echo "BLAS_LIBS   = " ${BLAS_LIBS} >> ${TARGET}
	echo "pthead_lib  = " ${pthread_lib} >> ${TARGET}
	echo "FLIBS       = " ${FLIBS} >> ${TARGET}
	echo "MUMPS_LIBS = "  ${MUMPS_LIBS} >> ${TARGET}
	echo "SDPA_LIBS = " ${SDPA_LIB} ${MUMPS_LIBS}\
	       ${LAPACK_LIBS} ${BLAS_LIBS} \
	       ${pthread_lib} ${FLIBS}  >> ${TARGET}

force: clean all
clean:
	rm -f make.inc *~
distclean: clean
install: all
uninstall:

