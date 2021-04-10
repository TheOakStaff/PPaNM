#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main(int argc, char const *argv[]) {
  if (argc<2) {
    printf("There was no arguments given to %s\n", argv[0]);
  }
  else {
    printf("read-cmdline\n");
    printf("d\tsin(d)\tcos(d)\n");
    for (int i = 1; i < argc; i++) {
      int d = atoi(argv[i]);
      printf("%d\t%.5f\t%.5f\n", d, sin(d), cos(d));
    }
  }

  return 0;
}
