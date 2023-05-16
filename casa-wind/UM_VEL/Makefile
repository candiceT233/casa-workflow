SOURCES.c= run_UM_VEL.c UM_VEL.c
$(CC)=gcc
CFLAGS=-Wall -g -I/usr/include/libxml2
LDFLAGS=-lnetcdf -lconfig -lm -lxml2

PROGRAM=UM_VEL
OBJECTS= $(SOURCES.c:.c=.o)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)


install: all
	@echo "nothing to install"
