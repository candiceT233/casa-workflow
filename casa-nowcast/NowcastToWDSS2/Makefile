#
## Makefile
#
SOURCES.c= run_convert.c NowcastToWDSS2.c
$(CC)=gcc 

CFLAGS=-Wall -g 
LDFLAGS=-lm -lnetcdf

PROGRAM=NowcastToWDSS2 

OBJECTS= $(SOURCES.c:.c=.o)

.KEEP_STATE:

debug := CFLAGS= -g -s

all debug: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)

