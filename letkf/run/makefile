#
# sbPOM makefile
#

#-----------------------------------------------------------------------
# Settings that depend on the system and the compiler
#-----------------------------------------------------------------------
# Set macros
CPP = cpp -P

### Fugaku ###
FC=mpifrtpx
LD=mpifrtpx

CLEAN = rm

# Set libraries and include files
#Fugaku
#NETCDFINC = -I/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/netcdf-c-4.7.3-t6k77xsvgeas4eaesj35f4ws6qy6ptwd/include
 
#NETCDFFINC = -I/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/netcdf-fortran-4.5.2-onfzr4kdmdychizfuaeyqul6ylaqpstv/include
#NETCDFLIB =  -L/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/netcdf-c-4.7.3-t6k77xsvgeas4eaesj35f4ws6qy6ptwd/lib
#NETCDFFLIB = -L/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/netcdf-fortran-4.5.2-onfzr4kdmdychizfuaeyqul6ylaqpstv/lib
#HDF5INC = -I/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/hdf5-1.10.6-r4ahikzbh42ezqqh5s7mntf2plinxt5s/include
#HDF5LIB = -L/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel8-a64fx/fj-4.2.0/hdf5-1.10.6-r4ahikzbh42ezqqh5s7mntf2plinxt5s/lib

#FFLAGS = -Kfast,SVE,parallel -Kopenmp -Nfjomplib -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast,SVE,parallel -Kopenmp -Nfjomplib -SCALAPACK -SSL2 $(NETCDFFLIB) -lnetcdff $(NETCDFLIB) -lnetcdf -lnetcdf -lm $(HDF5LIB) -lhdf5 -lhdf5_hl -lfjprofmpi

#JSS3
NETCDFINC = -I/opt/JX/oss/aarch64/netcdf/4.7.3/include
NETCDFFINC = -I/opt/JX/oss/aarch64/netcdf-fortran/4.5.2/include
NETCDFLIB = -L/opt/JX/oss/aarch64/netcdf/4.7.3/lib
NETCDFFLIB = -L/opt/JX/oss/aarch64/netcdf-fortran/4.5.2/lib 
HDF5INC = -I/opt/JX/oss/aarch64/hdf5/1.12.0/include
HDF5LIB = -L/opt/JX/oss/aarch64/hdf5/1.12.0/lib

FFLAGS = -Kfast,parallel,openmp -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = -Kfast,parallel,openmp -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) -lnetcdff ${HDF5LIB} $(NETCDFLIB) -lnetcdf -lnetcdf -lm -lhdf5_hl -lhdf5 -lm -lcurl

#FFLAGS = -Kfast -Kopenmp -Kparallel -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast -Kopenmp -Kparallel -SCALAPACK -SSL2 $(NETCDFFLIB) -lnetcdff ${HDF5LIB} $(NETCDFLIB) -lnetcdf -lnetcdf -lm -lhdf5_hl -lhdf5 -lm -lcurl


#-----------------------------------------------------------------------
# Set the executable
#-----------------------------------------------------------------------
BIN = letkf.exe

#-----------------------------------------------------------------------
# Define source directory
#-----------------------------------------------------------------------
SRCDIR = src

#-----------------------------------------------------------------------
# Define objects
#-----------------------------------------------------------------------
OBJS = SFMT.o \
       common_setting.o  \
       common.o \
       common_mpi.o \
       common_mtx.o \
       netlib.o \
       common_letkf.o  \
       common_pom.o \
       common_obs_pom.o  \
       common_mpi_pom.o  \
       letkf_obs.o \
       letkf_tools.o  \
       letkf.o

VPATH = $(SRCDIR)

#-----------------------------------------------------------------------
# Set implicit rules for compilation
#-----------------------------------------------------------------------
%.o: %.f
	@echo
	$(FC) -c $(FFLAGS) $<

%.o: %.f90
	@echo
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90 %.o
	@echo

#-----------------------------------------------------------------------
# Set implicit rules for dependencies
#-----------------------------------------------------------------------
%.f: %.F
	@echo
	$(CPP) $(FFLAGS) $< > $*.f

#-----------------------------------------------------------------------
# Create the executable
#-----------------------------------------------------------------------
$(BIN): $(OBJS)
	@echo
	$(LD) -o $(BIN) $(FFLAGS) $(OBJS) $(LIBS)

#-----------------------------------------------------------------------
# Cleaning target
#-----------------------------------------------------------------------
clean:
#	cd $(SRCDIR)
	@rm -f *.o *.mod
