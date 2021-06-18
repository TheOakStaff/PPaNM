#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"


void plainmc(
  int dim,
  double f(int dim, double *x),
  double *a,
  double *b,
  int N,
  gsl_vector* vec,
  int seed)
  {
        double V=1; for(int i=0;i<dim;i++){
          V*=b[i]-a[i];
        }
        double sum=0;
        double sum2=0;
        double x[dim];
        for(int i=0;i<N;i++){
          for(int i=0;i<dim;i++){
            x[i]=a[i]+ ((double)rand_r(&seed) / (double)RAND_MAX)*(b[i]-a[i]);
          }
          double fx=f(dim,x);
          sum+=fx;
          sum2+=fx*fx;
          }
        double mean=sum/N;
        double sigma=sqrt(sum2/N-mean*mean);
        gsl_vector_set(vec,0,mean*V);
        gsl_vector_set(vec,1,sigma*V/sqrt(N));
}

int main(int argc, char const *argv[]) {
  int dim = 3;
  double f(int dim, double *x){
    return 1/((1 - cos(x[0]) * cos(x[1]) * cos(x[2])) * M_PI * M_PI * M_PI);
  }
  double a[] = {0,0,0};
  double b[] = {M_PI,M_PI,M_PI};
  int N = 1e6;
  gsl_vector* vec = gsl_vector_alloc(dim);
  int seed = 5;

  plainmc(dim,f,a,b,N,vec,seed);
  vector_print(stdout,vec);
  printf("∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]⁻¹ = 1.3932039296856768591842462603255\n");
  printf("Monte Carlo integration with 1e6 points retuns %.10g\n",gsl_vector_get(vec,0));
  printf("with error = %.10g\n\n", gsl_vector_get(vec,1));

  int dim2 = 3;
  double f2(int dim2, double *x2){
    return x2[0]*x2[0] + 2*x2[1]*x2[1] + 3*x2[2]*x2[2];
  }
  double a2[] = {0,0,0};
  double b2[] = {3,2,1};
  int N2 = 1e7;
  gsl_vector* vec2 = gsl_vector_alloc(dim2);
  int seed2 = 2;

  plainmc(dim2,f2,a2,b2,N2,vec2,seed2);
  vector_print(stdout,vec2);

  printf("∫03  dx ∫02  dy ∫01  dz x²+2y²+3z² = 40\n");
  printf("Monte Carlo integration with 1e7 points retuns %.10g\n",gsl_vector_get(vec2,0));
  printf("with error = %.10g\n\n", gsl_vector_get(vec2,1));

  // PART B

  return 0;
}
