#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main() {
  printf("\n\n\n");
  printf("Equal test function\n");

  int equal(double a, double b, double tau, double epsilon) {
    printf("running Equal test with a = %g, b = %g, tau = %g and epsilon = % g\n",a,b,tau,epsilon);
    if (sqrt((a - b) * (a - b)) < tau ) {
      printf("Absolute pressision achived\n");
      return 1;
    }
    else if (sqrt((a - b) * (a - b)) / (sqrt(a * a) + sqrt(b * b)) < epsilon / 2) {
      printf("Relative pressision achived\n");
      return 1;
    }
    else {
      printf("Numbers are not equal\n");
      return 0;
    }
  };

  equal(2.000,2.002,0.001,0.001);

  return 0;
}
