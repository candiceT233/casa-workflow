SOURCES.c= run_mvt.c mvt.c threshold_functions.c
$(CC)=gcc
CFLAGS=-Wall -I. -I/usr/local/include -I/usr/include/libxml2 -I/usr/local/jansson/include -I/usr/local/libconfig/include
LDLIBS=-L/usr/local/jansson/lib -L/usr/local/libconfig/lib -lm -lnetcdf -ljansson -lconfig -lxml2

PROGRAM=mvt
OBJECTS= $(SOURCES.c:.c=.o)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:
	rm -f $(PROGRAM) $(OBJECTS)

install: all
	@echo "nothing to install"
