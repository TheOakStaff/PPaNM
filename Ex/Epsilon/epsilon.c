#include "stdio.h"
#include "math.h"
#include "limits.h"
#include "float.h"

int main() {

  int i=1; while (i+1>i) {
    i++;
  }
  printf("my max int = %i\n", i);

  i = 1;
  for (i = 1; i < i+1; i++) {
  }
  printf("my max int = %i\n", i);

  i = 1;
  do {
    i++;
  } while(i+1>1);
  printf("my max int = %i\n", i);

  return 0;
}
