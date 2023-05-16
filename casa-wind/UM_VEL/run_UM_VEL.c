#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <argp.h>
#include <math.h>
#include <sys/resource.h>

char commandline[]= "UM_VEL [input_filename0] ... [input_filename_n]";

int main(int argc, char* argv[])
{

  double t1,t2;
  t1=clock();

  if(argc == 1) {
    printf("Usage: %s\n", commandline);
    exit(-1);
  }
  
  const rlim_t sSz = 64L * 1024L * 1024L;   // min stack size = 64 Mb
  struct rlimit rlim;
  int stat;

  stat = getrlimit(RLIMIT_STACK, &rlim);
  if (stat == 0) {
    if (rlim.rlim_cur < sSz) {
      rlim.rlim_cur = sSz;
      stat = setrlimit(RLIMIT_STACK, &rlim);
      if (stat != 0) {
	fprintf(stderr, "setrlimit returned result = %d\n", stat);
      }
    }
  }

  /*Function Declaration */
  void UM_VEL(size_t count, char* i_filenames[]);
  
  /* Function call */
  UM_VEL(argc - 1, argv + 1 );
  
  t2=clock();

  return 0;
}






