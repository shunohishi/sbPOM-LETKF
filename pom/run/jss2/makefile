#
# sbPOM makefile
#

#-----------------------------------------------------------------------
# Settings that depend on the system and the compiler
#-----------------------------------------------------------------------
# Set macros
CPP = cpp -P

### JSS2 ###
FC=f90sjx
LD=f90sjx

CLEAN = rm

# Set libraries and include files
### JSS2 ###
NETCDFINC = -I/opt/JX/oss/netcdf-4.2.1.1/s64fx/include
NETCDFFINC = -I/opt/JX/oss/netcdf-fortran-4.2/s64fx/include
NETCDFLIB = -L/opt/JX/oss/netcdf-4.2.1.1/s64fx/lib
NETCDFFLIB = -L/opt/JX/oss/netcdf-fortran-4.2/s64fx/lib
HDF5INC = -I/opt/JX/oss/hdf5-1.8.9/s64fx/include
HDF5LIB = -L/opt/JX/oss/hdf5-1.8.9/s64fx/lib
#CHECK = -Haefosux

### JSS2 ###
FFLAGS = -Umpi -O3 $(CHECK) $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = $(NETCDFFLIB) -lnetcdff $(NETCDFLIB) -lnetcdf $(HDF5LIB) -lhdf5_hl -lhdf5 -lz

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
