#
## Makefile
#
SOURCES.c= pointAlert.c
$(CC)=gcc 
CFLAGS=-Wall -g -I/usr/local/include -I/usr/local/jansson/include -I/usr/local/libconfig/include
LDLIBS= -L/usr/local/jansson/lib -L/usr/local/libconfig/lib -lm -ljansson -lconfig

PROGRAM=pointAlert

OBJECTS= $(SOURCES.c:.c=.o)

.KEEP_STATE:

debug := CFLAGS= -g -s

all debug: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)

