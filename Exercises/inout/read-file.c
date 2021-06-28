#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main(int argc, char const *argv[]) {
  FILE* my_in_stream = fopen(argv[1],"r");
  FILE* my_out_stream = fopen(argv[2],"w");
  int items = 1;
  int x;
  fprintf(my_out_stream,"read-file\n");
  fprintf(my_out_stream,"d\tsin(d)\tcos(d)\n");
  while (items == 1) {
    items = fscanf(my_in_stream,"%d\n",&x);
    fprintf(my_out_stream,"%d\t%.5f\t%.5f\n", x, sin(x), cos(x));
  }
  fclose(my_in_stream);
  fclose(my_out_stream);
  return 0;
}
