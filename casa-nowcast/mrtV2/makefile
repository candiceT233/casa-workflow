SOURCES.c= run_mrtV2.c mrtV2.c threshold_functions.c

$(CC)=gcc
# CFLAGS=-Wall -I. -I/usr/local/include -I/usr/include/libxml2 -I/usr/local/libconfig/include -I/usr/local/jansson/include
# LDLIBS= -L/usr/local/libconfig/lib -L/usr/local/jansson/lib -lm  -lnetcdf -ljansson -lconfig -lxml2


# LIBCONFIG_PATH := $(shell spack find -p libconfig | grep -oE '/[^[:space:]]+')
# JANSSON_PATH := $(shell spack find -p jansson | grep -oE '/[^[:space:]]+')
# LIBXML2_PATH := $(shell spack find -p libxml2 | grep -oE '/[^[:space:]]+')
# LIBXML2_INCLUDE=/usr/include/libxml2

LIBCONFIG_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libconfig-1.7.2-htzxtlcxpvlcyb3cti7sr7ttnzuiwdjb
JANSSON_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/jansson-2.9-6mxyawnx3ftbm325dyjlumi3m2bw2c3a
LIBXML2_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libxml2-2.10.1-vfgoey2prebiqr3lljtvdoun4sciczsf

NETCDF_PATH=/qfs/people/tang584/install/netcdf-c
HDF5_PATH=/qfs/people/tang584/install/hdf5

SYSTEM_INCLUDE=/usr/include
SYSMTE_LIB64=/usr/lib64
SYSTEM_LIB=/usr/lib

CFLAGS=-Wall -I. -I$(LIBXML2_PATH)/include -I$(LIBXML2_PATH)/include/libxml2 -I$(NETCDF_PATH)/include -I$(HDF5_PATH)/include -I$(LIBCONFIG_PATH)/include -I$(JANSSON_PATH)/include -I$(SYSTEM_INCLUDE) 
LDLIBS=-L$(LIBCONFIG_PATH)/lib -L$(LIBXML2_PATH)/lib -L$(JANSSON_PATH)/lib -L$(NETCDF_PATH)/lib -L$(HDF5_PATH)/lib #-L$(SYSMTE_LIB64) #-L$(SYSTEM_LIB)

LDFLAGS=-Wl,-rpath,$(LIBCONFIG_PATH)/lib $(LDLIBS) -lm -lnetcdf -lhdf5 -lhdf5_hl -ljansson -lconfig -lxml2 -lcurl

PROGRAM=mrtV2
OBJECTS= $(SOURCES.c:.c=.o)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
	rm -f $(PROGRAM) $(OBJECTS)

checkp:
	@echo "LIBCONFIG_PATH: $(LIBCONFIG_PATH)"
	@echo "LIBXML2_PATH: $(LIBXML2_PATH)"
	@echo "JANSSON_PATH: $(JANSSON_PATH)"
	@echo "NETCDF_PATH: $(NETCDF_PATH)"
	@echo "HDF5_PATH: $(HDF5_PATH)"

install: all
	@echo "nothing to install"
