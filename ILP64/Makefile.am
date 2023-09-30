SUBDIRS = etc mumps

exec_prefix    = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/
includedir     = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/include
datadir        = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/share
etcdir         = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/etc
docdir         = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/doc
infodir        = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/info
localstatedir  = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/var
mandir         = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/man
sharedstatedir = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/share
sysconfdir     = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/etc
exampledir     = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/example
libexampledir  = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/libexample
mexdir         = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/mex

# each file of mumps is copied by 'nobase' option
#   due to sub-directory
mumpsdir     = ${prefix}/@PACKAGE_NAME@.@PACKAGE_VERSION@/


etc_DATA = etc/make.inc
doc_DATA = doc/sdpa7-manual.pdf

include_HEADERS = sdpa_algebra.h sdpa_block.h \
sdpa_call.h sdpa_chordal.h sdpa_dataset.h sdpa_dpotrf.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h

lib_LIBRARIES = libsdpa.a
libsdpa_a_SOURCES = sdpa_block.cpp sdpa_call.cpp sdpa_chordal.cpp \
sdpa_dataset.cpp sdpa_dpotrf.cpp sdpa_io.cpp sdpa_jordan.cpp \
sdpa_linear.cpp sdpa_newton.cpp sdpa_parts.cpp sdpa_solve.cpp \
sdpa_struct.cpp sdpa_tool.cpp \
sdpa_algebra.h sdpa_block.h \
sdpa_call.h sdpa_chordal.h sdpa_dataset.h sdpa_dpotrf.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h  
libsdpa_a_CXXFLAGS = $(pthread_cflags) -fPIC -I$(MUMPS_DIR)/include

bin_PROGRAMS = sdpa
sdpa_SOURCES = sdpa_exe.cpp
sdpa_LDADD = $(pthread_lib) -L. -lsdpa -L$(MUMPS_DIR)/lib  -ldmumps -lmumps_common -lpord -L$(MUMPS_DIR)/libseq -lmpiseq $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 
sdpa_CXXFLAGS = $(pthread_cflags)  -fPIC -I$(MUMPS_DIR)/include

nobase_mumps_HEADERS = mumps/build/include/dmumps_c.h \
mumps/build/include/dmumps_root.h \
mumps/build/include/dmumps_struc.h \
mumps/build/include/mumps_c_types.h \
mumps/build/include/mumps_compat.h

nobase_mumps_DATA = mumps/build/lib/libdmumps.a \
mumps/build/lib/libmumps_common.a \
mumps/build/lib/libpord.a \
mumps/build/libseq/libmpiseq.a \
bin/param.sdpa
# param.sdpa should installed into bin/

example_DATA = example1.dat-s example1.dat example1.ini \
example1.ini-s example2.dat

libexample_DATA  = libexample/Makefile \
libexample/example1.cpp libexample/example2.cpp \
libexample/example3.cpp libexample/example4.cpp \
libexample/example5.cpp libexample/example6.cpp

mex_DATA = mex/Makefile mex/gensdpafile.m mex/initial_point.m \
mex/mexSedumiWrap.cpp mex/mexsdpa.cpp \
mex/mexAggSDPcones.cpp mex/mexDisAggSDPsol.cpp mex/mexFprintf.c \
mex/param.m mex/read_data.m mex/sdpam.m mex/sedumiwrap.m 


