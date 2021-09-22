
#######
# shell
#######
SHELL = /bin/sh

F90 = gfortran

BASEDIR		=	/Users/joshua/bob/
HDF5BASEDIR	=	$(BASEDIR)/usr/hdf5/
MPIBASEDIR	=	$(BASEDIR)/usr/openmpi/

#wow wtf caps matters below (shell not SHELL)
#this may not work					\here/ it does, just with 'shell'
MPI_FCOMPILEFLAGS	=	$(shell /Users/joshua/bob/usr/openmpi/bin/mpif77 --showme:compile) \
						$(shell /Users/joshua/bob/usr/openmpi/bin/mpif90 --showme:compile)
MPI_FLINKFLAGS		=	$(shell $(MPIBASEDIR)/bin/mpif77 --showme:link) \
						$(shell $(MPIBASEDIR)/bin/mpif90 --showme:link)

H5_FCOMPILEFLAGS = -I$(HDF5BASEDIR)/lib -I$(HDF5BASEDIR)/include
H5_FLINKFLAGS = -L$(HDF5BASEDIR)/lib -lhdf5_fortran -lhdf5


COMPILATION_TYPE = production #config/osiris_config

#I haven't checked these for this system, probably just leave out
F90FLAGS_production = -cpp -O3 -ipo -no-prec-div


F90F = $(F90FLAGS_$(COMPILATION_TYPE)) $(MPI_FCOMPILEFLAGS) $(H5_FCOMPILEFLAGS) \
	$(H5_FLINKFLAGS) $(MPI_FLINKFLAGS)

#F90F = $(MPI_FCOMPILEFLAGS) $(MPI_FLINKFLAGS)


h5_util-nohash.o	:	h5_util-nohash.f90
	$(F90) $(F90F) -c h5_util-nohash.f90