#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/* --------------------------------------------------------------------------*/
void stat_dir(char directory[2048], int *n) {

  struct stat sb;

  /*printf("%s \n",directory);*/
  
  if (stat(directory, &sb) == -1) {
    *n = -1;
  } else {
    *n = 0;
  }
  
}

/* --------------------------------------------------------------------------*/
void make_dir(char directory[2048], int *n) {

  /*printf("%s \n",directory);*/

  if (mkdir(directory, S_IRWXU ) == -1) {
      *n = -1;
  } else {
    *n = 0;
  }
  
}
