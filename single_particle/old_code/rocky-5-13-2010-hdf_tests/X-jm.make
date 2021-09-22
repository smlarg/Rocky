#Makefile for 1D PIC codes in new_beps1.source

GOBJS = nullgks1.o
CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile Intel compiler with Mac OS X

FC90 = gfortran
FC77 = gfortran
CC = gcc

OPTS90 = -O3 -r8
OPTS77 = -O3 -r8
CCOPTS = -O
MOPTS = -save
MBOPTS = -save
LOPTS = -lpthread
NOSAVE = -automatic
LEGACY =

HDF_DIR =  /Users/joshua/bob/
H5_DIR =   /Users/joshua/bob/
SZ_DIR =   /Users/joshua/bob/
JPEG_DIR = /Users/joshua/bob/

INCPATH = -I$(HDF_DIR)/include -I$(H5_DIR)/include -I$(H5_DIR)/lib 

LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
	-L$(HDF_DIR)/lib -lmfhdf -ldf -ljpeg -lz \
	-L$(H5_DIR)/lib -lsz -lhdf5_fortran -lhdf5 \
	-L$(SZ_DIR)/lib \
	-L$(JPEG_DIR)/lib


ESOBJS = globals.o init1mod.o diag1mod.o bpush1mod.o push1mod.o fft1mod.o \
field1mod.o init1lib.o bpush1lib.o push1lib.o fft1lib.o field1lib.o diag1lib.o \
ext_driver_jf.o init1mod_jf.o diag_jf.o ampere_jf.o hdf_write_nompi_jf.o \
bpush1mod.o bpush1lib.o

EMOBJS = dpush1mod.o dpush1lib.o

# Linkage rule

all : h5ex_d_compact.e

#all : new_beps1_jf.out
#all : new_beps1.out new_bbeps1.out new_dbeps1.out


new_beps1.out : new_beps1.o $(ESOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps1.out \
        new_beps1.o $(ESOBJS) $(GOBJS) $(LIBS)

h5ex_d_compact.e : h5ex_d_compact.o $(ESOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) mkdir_f_fxns.o -o new_beps1_jf.out \
  h5ex_d_compact.o $(ESOBJS) $(GOBJS) $(LIBS) $(INCPATH)

#new_beps1_jf.out : new_beps1_jf.o $(ESOBJS) $(GOBJS)
#	$(FC90) $(OPTS90) $(LOPTS) mkdir_f_fxns.o -o new_beps1_jf.out \
#  new_beps1_jf.o $(ESOBJS) $(GOBJS) $(LIBS) $(INCPATH)



new_bbeps1.out : new_bbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_bbeps1.out \
        new_bbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS) $(LIBS)

new_dbeps1.out : new_dbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_dbeps1.out \
        new_dbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS) $(LIBS)

new_beps1gl.out : new_beps1gl.o $(ESOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps1gl.out \
        new_beps1gl.o $(ESOBJS) $(GOBJS) $(LIBS)

# Compilation rules

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libpsp.o : libpsp.f
	$(FC77) $(OPTS77) -c libpsp.f

librstr.o : librstr.f
	$(FC77) $(OPTS77) -c librstr.f

ncarstub.o : ncarstub.f
	$(FC77) $(OPTS77) -c ncarstub.f

libplt10.o : libplt10.f
	$(FC77) $(OPTS77) -c libplt10.f

plot10.o : plot10.f
	$(FC77) $(OPTS77) -c plot10.f

libt1.o : libt1.f
	$(FC77) $(OPTS77) -c libt1.f

libloc1.o : libloc1.f
	$(FC77) $(OPTS77) -c libloc1.f

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f

init1lib.o : init1lib.f
	$(FC77) $(OPTS77) -c init1lib.f

push1lib.o : push1lib.f
	$(FC77) $(OPTS77) -c push1lib.f

bpush1lib.o : bpush1lib.f
	$(FC77) $(OPTS77) -c bpush1lib.f

dpush1lib.o : dpush1lib.f
	$(FC77) $(OPTS77) -c dpush1lib.f

fft1lib.o : fft1lib.f
	$(FC77) $(OPTS77) -c fft1lib.f

field1lib.o : field1lib.f
	$(FC77) $(OPTS77) -c field1lib.f

diag1lib.o : diag1lib.f
	$(FC77) $(OPTS77) -c diag1lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

init1mod.o : init1mod.f globals.o
	$(FC90) $(OPTS90) -c init1mod.f

push1mod.o : push1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c push1mod.f

bpush1mod.o : bpush1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c bpush1mod.f

dpush1mod.o : dpush1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c dpush1mod.f

fft1mod.o : fft1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c fft1mod.f

field1mod.o : field1mod.f globals.o
	$(FC90) $(OPTS90) $(LEGACY) -c field1mod.f

diag1mod.o : diag1mod.f init1mod.o
	$(FC90) $(OPTS90) -c diag1mod.f

mkdir_f_fxns.o : mkdir_f_fxns.c
	$(CC) $(CCOPTS) -c mkdir_f_fxns.c

init1mod_jf.o : init1mod_jf.f
	$(FC90) $(OPTS90) -free -c init1mod_jf.f

ext_driver_jf.o : ext_driver_jf.f init1mod_jf.o
	$(FC90) $(OPTS90) -free -c ext_driver_jf.f

diag_jf.o : diag_jf.f init1mod.o init1mod_jf.o hdf_write_nompi_jf.o
	$(FC90) $(OPTS90) $(INCPATH) -free -c diag_jf.f
	
ampere_jf.o : ampere_jf.f
	$(FC90) $(OPTS90) -free -c ampere_jf.f

hdf_write_nompi_jf.o : hdf_write_nompi_jf.f
	$(FC90) $(OPTS90) $(INCPATH) -c -free hdf_write_nompi_jf.f

new_beps1.o : new_beps1.f push1mod.o fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps1.f

new_beps1_jf.o : new_beps1_jf.f push1mod.o fft1mod.o field1mod.o init1mod_jf.o \
              ext_driver_jf.o diag_jf.o ampere_jf.o hdf_write_nompi_jf.o \
              mkdir_f_fxns.o bpush1mod.o
	$(FC90) $(OPTS90) $(MOPTS) $(INCPATH) -free -c new_beps1_jf.f

new_dbeps1.o : new_dbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
               fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_dbeps1.f

new_bbeps1.o : new_bbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
               fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_bbeps1.f

new_beps1gl.o : new_beps1gl.f push1mod.o fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps1gl.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f *.out
