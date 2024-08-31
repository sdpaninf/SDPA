MUMPS_VER = 5.7.0

SDPA_MAKEFILE = ../Makefile
# This file is automatically made by './configure' script
###################################################################
# The following variables are set from
# the arguments given to './configure' script for building SDPA
###################################################################

CC          = ${shell grep "^CC ="        ${SDPA_MAKEFILE} | cut -f2 -d=}
CFLAGS      = ${shell grep "^CFLAGS = "   ${SDPA_MAKEFILE} | cut -f2 -d=}
CXX         = ${shell grep "^CXX ="       ${SDPA_MAKEFILE} | cut -f2 -d=}
CXXFLAGS    = ${shell grep "^CXXFLAGS = " ${SDPA_MAKEFILE} | cut -f2 -d=}
F77         = ${shell grep "^F77 ="       ${SDPA_MAKEFILE} | cut -f2 -d=}
FFLAGS      = ${shell grep "^FFLAGS = "   ${SDPA_MAKEFILE} | cut -f2 -d=}
###################################################################

all: build/lib/libdmumps.a

build/lib/libdmumps.a:
	rm -rf build;
	tar xzf archive/MUMPS_${MUMPS_VER}.tar.gz
	mv MUMPS_${MUMPS_VER} build;
	cd build; cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc;

	cd build; echo "CC = "   ${CC}       >> Makefile.inc;
	cd build; echo "FC = "   ${F77}      >> Makefile.inc;
	cd build; echo "FL = "   ${F77}      >> Makefile.inc;
	cd build; echo "OPTC = " ${CFLAGS}   >> Makefile.inc;
	cd build; echo "OPTF = " ${FFLAGS}   >> Makefile.inc;
	cd build; echo "OPTL = " ${FFLAGS}   >> Makefile.inc;
	cd build; echo "LPORDDIR = "`pwd`"/PORD/lib/"     >> Makefile.inc;
	cd build; echo "IPORDDIR = "`pwd`"/PORD/include/" >> Makefile.inc;
	cd build ; ${MAKE} d

clean:
	@echo "##########################################################"
	@echo "If you want to rebuild mumps dir, use 'make distclean'"
	@echo "##########################################################"

distclean: 
	rm -rf build

install: all
uninstall:

