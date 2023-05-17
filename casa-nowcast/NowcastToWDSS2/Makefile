#
## Makefile
#
SOURCES.c= run_convert.c NowcastToWDSS2.c
$(CC)=gcc 

LIBXML2_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libxml2-2.10.1-vfgoey2prebiqr3lljtvdoun4sciczsf
NETCDF_PATH=/qfs/people/tang584/install/netcdf-c
HDF5_PATH=/qfs/people/tang584/install/hdf5

LDLIBS=-L$(NETCDF_PATH)/lib -L$(HDF5_PATH)/lib -L$(LIBXML2_PATH)/lib

CFLAGS=-Wall -g -I$(HDF5_PATH)/include -I$(NETCDF_PATH)/include -I$(LIBXML2_PATH)/include/libxml2
LDFLAGS=$(LDLIBS) -lm -lhdf5 -lhdf5_hl -lnetcdf -lxml2 -lcurl

PROGRAM=NowcastToWDSS2 

OBJECTS= $(SOURCES.c:.c=.o)

.KEEP_STATE:

debug := CFLAGS= -g -s

all debug: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)

