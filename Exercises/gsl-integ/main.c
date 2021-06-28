#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>


double fun (double x, void * params) {
  // Create dummy variable. This is neede cause gsl is stupid
  double k = *(double*) params;
  double g = log(k*x) / sqrt(x);
  return g;
}

int main() {
  // First allocate memory as such. This is needed cause gsl is stupid.
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  // gsl is stupid and needs a double result and a double error
    double result, error;
    double k = 1.0;

    // This is neede cause... gsl says so.
    gsl_function F;
    F.function = &fun;
    F.params = &k;

    gsl_integration_qags(&F,0,1,0,1e-6,1000,w,&result,&error);

    printf("result = %.9f\n",result);

  return 0;
}
