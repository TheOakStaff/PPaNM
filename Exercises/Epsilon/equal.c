#include <stdio.h>
#include <math.h>

double sqrt(double x);

int equal(double a, double b, double tau, double epsilon) {
  printf("\n\n\n");
  printf("Equal test function\n");
  printf("running Equal test with a = %g, b = %g, tau = %g and epsilon = % g\n",a,b,tau,epsilon);
  double absolute = fabsf(a - b);
  double relative = fabsf(a - b) / (fabsf(a) + fabsf(b));
  if ( absolute< tau ) {
    printf("Absolute pressision achived\n");
    return 1;
  }
  else if (relative < epsilon / 2) {
    printf("Relative pressision achived\n");
    return 1;
  }
    else {
      printf("Numbers are not equal\n");
    return 0;
  }
};
