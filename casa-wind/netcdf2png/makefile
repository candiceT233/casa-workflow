
# CFLAGS=-Wall -O -I/usr/local/netcdf/include -I/usr/include/libxml2
# LDLIBS=-lm -lxml2 -lz -lpng -L/usr/local/netcdf/lib -lnetcdf 

CC=gcc

LIBXML2_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libxml2-2.10.1-vfgoey2prebiqr3lljtvdoun4sciczsf
NETCDF_PATH=/qfs/people/tang584/install/netcdf-c
HDF5_PATH=/qfs/people/tang584/install/hdf5

CFLAGS=-Wall -O -I$(HDF5_PATH)/include -I$(NETCDF_PATH)/include -I$(LIBXML2_PATH)/include/libxml2
LDLIBS=-lm -lz -lpng -L$(NETCDF_PATH)/lib -L$(HDF5_PATH)/lib -L$(LIBXML2_PATH)/lib -lnetcdf -lhdf5 -lhdf5_hl -lxml2 -lcurl


all: merged_netcdf2png

install: all
	@echo "nothing to install"
