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

#FFLAGS = -Kfast,parallel -Kopenmp -Nfjomplib $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast,parallel -Kopenmp -Nfjomplib $(NETCDFFLIB) -lnetcdff $(NETCDFLIB) -lnetcdf -lnetcdf -lm $(HDF5LIB) -lhdf5 -lhdf5_hl -lfjprofmpi

#JSS3
NETCDFINC = -I/opt/JX/oss/aarch64/netcdf/4.7.3/include
NETCDFFINC = -I/opt/JX/oss/aarch64/netcdf-fortran/4.5.2/include
NETCDFLIB = -L/opt/JX/oss/aarch64/netcdf/4.7.3/lib
NETCDFFLIB = -L/opt/JX/oss/aarch64/netcdf-fortran/4.5.2/lib 
HDF5INC = -I/opt/JX/oss/aarch64/hdf5/1.12.0/include
HDF5LIB = -L/opt/JX/oss/aarch64/hdf5/1.12.0/lib

FFLAGS = -Kfast -Kopenmp -Kparallel -Nalloc_assign $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = -Kfast -Kopenmp -Kparallel -Nalloc_assign $(NETCDFFLIB) -lnetcdff ${HDF5LIB} $(NETCDFLIB) -lnetcdf -lnetcdf -lm -lhdf5_hl -lhdf5 -lm -lcurl

#-----------------------------------------------------------------------
# Set the executable
#-----------------------------------------------------------------------
BIN = pom.exe

#-----------------------------------------------------------------------
# Define source directory
#-----------------------------------------------------------------------
SRCDIR = pom

#-----------------------------------------------------------------------
# Define objects
#-----------------------------------------------------------------------
OBJS = mod_const.o  \
       mod_thermo.o \
       mod_blk_coare.o \
       mod_blk_ecmwf.o \
       mod_blk_ncar.o  \
       mod_aerobulk_compute.o \
       mod_aerobulk.o  \
       parallel_mpi.o  \
       io_netcdf_relax.o \
       merge_scatter.o  \
       timelib.o        \
       daily_average.o  \
       fluxlib_rh_rt.o  \
       solver.o         \
       turb_MYNNF_lev2.5_v20130905.o \
       profq_mynnf.o    \
       bounds_forcing.o \
       bounds_forcing_frs.o \
       advance_relax.o      \
       initialize_relax.o   \
       iomngr_v20090918.o   \
       utils_iomngr_20090129.o \
       caldat_m.o         \
       com_serv_v201002.o \
       tide_serv_v20110624.o     \
       time_serv.o               \
       tide_jcope_v20110915.o    \
       interp2d_bispline_all_v201001.o   \
       com_mask_correction_all_v201001.o \
       tide_sbpom.o              \
       pom.o
#       uvstore.o                 \
#       uvstatistics.o            \

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
	@rm -f *.o *.mod *~ *exe
