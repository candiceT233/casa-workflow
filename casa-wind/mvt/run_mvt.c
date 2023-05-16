#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <argp.h>
#include "threshold_functions.h"

char commandline[]= "mvt [input_filename]" ;

int main(int argc, char* argv[])
{
  char *i_filename;

  if(argc != 2) {
    printf("Usage: %s\n", commandline);
    exit(-1);
  }

  i_filename = argv[1];

  mvt(i_filename);

  return 0;
}
