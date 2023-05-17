#
## Makefile
#
SOURCES.c= pointAlert.c
$(CC)=gcc 
# CFLAGS=-Wall -g -I/usr/local/include -I/usr/local/jansson/include -I/usr/local/libconfig/include
# LDLIBS= -L/usr/local/jansson/lib -L/usr/local/libconfig/lib -lm -ljansson -lconfig


LIBCONFIG_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libconfig-1.7.2-htzxtlcxpvlcyb3cti7sr7ttnzuiwdjb
JANSSON_PATH=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/jansson-2.9-6mxyawnx3ftbm325dyjlumi3m2bw2c3a

CFLAGS=-Wall -g -I$(LIBCONFIG_PATH)/include -I$(JANSSON_PATH)/include
LDLIBS=-L$(LIBCONFIG_PATH)/lib -L$(JANSSON_PATH)/lib -lm -ljansson -lconfig

LDFLAGS=-Wl,-rpath,$(LIBCONFIG_PATH)/lib $(LDLIBS)

PROGRAM=pointAlert

OBJECTS= $(SOURCES.c:.c=.o)

.KEEP_STATE:

debug := CFLAGS= -g -s

all debug: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)

