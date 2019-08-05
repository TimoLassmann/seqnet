# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@


ACLOCAL_AMFLAGS = -I m4




AM_CFLAGS =  @SIMD_FLAGS@ @OPENMP_CFLAGS@ 
CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa test_dna.lfa libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5 test1.h5 test2.h5

# top_builddir = @top_builddir@


AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef


MYLIBDIRS = ../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lpthread -lm 


#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)

AM_CPPFLAGS=${MYINCDIRS} ${HDF5_CPPFLAGS}



#LOGO_SOURCES = motif_logo.c motif_logo.h

bin_PROGRAMS = seqnet

seqnet_SOURCES = \
run_seqnet.c \
alphabet.h \
alphabet.c \
bpm.h \
bpm.c \
parameters.h \
parameters.c \
rwalign.c \
global.h \
msa.h \
bisectingKmeans.h \
bisectingKmeans.c \
euclidean_dist.h \
euclidean_dist.c \
pick_anchor.h \
pick_anchor.c \
sequence_distance.h \
sequence_distance.c \
matrix_io.h \
matrix_io.c



check_PROGRAMS =  bpm_test rwaln alphabet
TESTS = bpm_test
TESTS_ENVIRONMENT = $(VALGRIND)

rwaln_SOURCES = \
rwalign.c \
alphabet.h \
alphabet.c \
msa.h
rwaln_CPPFLAGS = $(AM_CPPFLAGS) -DRWALIGN_TEST

bpm_test_SOURCES = \
bpm.h \
bpm.c \
alphabet.h \
alphabet.c
bpm_test_CPPFLAGS = $(AM_CPPFLAGS) -DBPM_UTEST


alphabet_SOURCES = \
alphabet.h \
alphabet.c
alphabet_CPPFLAGS = $(AM_CPPFLAGS) -DITEST

