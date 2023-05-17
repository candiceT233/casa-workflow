SOURCES.c= run_UM_VEL.c UM_VEL.c
$(CC)=gcc
# CFLAGS=-Wall -g -I/usr/include/libxml2
# LDFLAGS=-lnetcdf -lconfig -lm -lxml2

LIBCONFIG_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libconfig-1.7.2-htzxtlcxpvlcyb3cti7sr7ttnzuiwdjb
LIBXML2_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libxml2-2.10.1-vfgoey2prebiqr3lljtvdoun4sciczsf
NETCDF_PATH=/qfs/people/tang584/install/netcdf-c
HDF5_PATH=/qfs/people/tang584/install/hdf5

CFLAGS=-Wall -g -I$(HDF5_PATH)/include -I$(NETCDF_PATH)/include -I$(LIBXML2_PATH)/include/libxml2 -I$(LIBCONFIG_PATH)/include
LDLIBS=-L$(LIBCONFIG_PATH)/lib -L$(NETCDF_PATH)/lib -L$(HDF5_PATH)/lib -L$(LIBXML2_PATH)/lib
LDFLAGS=-Wl,-rpath,$(LIBCONFIG_PATH)/lib $(LDLIBS) -lnetcdf -lhdf5 -lhdf5_hl -lconfig -lm -lxml2 -lcurl

PROGRAM=UM_VEL
OBJECTS= $(SOURCES.c:.c=.o)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)


install: all
	@echo "nothing to install"
