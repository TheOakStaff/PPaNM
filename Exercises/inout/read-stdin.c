#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main() {
  int items = 1;
  int x;
  printf("read-stdin\n");
  printf("d\tsin(d)\tcos(d)\n");
  while (items == 1) {
    items = scanf("%d\n",&x);
    printf("%d\t%.5f\t%.5f\n", x, sin(x), cos(x));
  }

  return 0;
}
