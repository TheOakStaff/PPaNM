#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"
#include <assert.h>


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

double corput(int id, int base){
  double corput_q = 0;
  double coprime_b = (double)1/base;
  while(id > 0){
    corput_q += (id%base)*coprime_b;
    id /= base;
    coprime_b /= base;
}
  return corput_q;
}

void halton(int id, int dim, double* p){
  int b[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79}; // Series of prime numbers
  int b_Dim = sizeof(b)/sizeof(int);
  assert(dim <= b_Dim );
  for(int i = 0; i < dim; i++){
    p[i] = corput(id + 1, b[i]);
  }
}

void halton2(int id, int dim, double* p){
  int b[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83}; // Series of prime numbers
  int b_Dim = sizeof(b)/sizeof(int);
  assert(dim <= b_Dim );
  for(int i = 0; i < dim; i++){
    p[i] = corput(id + 1, b[i]);
  }
}

void HCRandomP(int id, int dim, double* lowerBound, double* upperBound, double* pts){
  halton(id, dim, pts);
  for(int i = 0; i < dim; i++){
    pts[i] = lowerBound[i] + pts[i]*(upperBound[i] - lowerBound[i]);
  }
}

void HCRandomP2(int id, int dim, double* lowerBound, double* upperBound, double* pts){
  halton2(id, dim, pts);
  for(int i = 0; i < dim; i++){
    pts[i] = lowerBound[i] + pts[i]*(upperBound[i] - lowerBound[i]);
  }
}

void HC_MC(
  int dim,
  double f(int dim, double *x),
  double *a,
  double *b,
  int N,
  gsl_vector* vec)
  {
        double V=1; for(int i=0;i<dim;i++){
          V*=b[i]-a[i];
        }
        double sum=0;
        double sum2=0;
        double sumS=0;
        double sum2S=0;
        double x[dim];
        double x2[dim];
        for(int i=0;i<N/2;i++){
          HCRandomP(i,dim,a,b,x);
          HCRandomP2(i,dim,a,b,x2);
          double fx=f(dim,x);
          double fx2=f(dim,x2);
          if(!isinf(fx) && !isinf(fx2)){
            sum +=fx;
            sumS +=fx*fx;
            sum2 +=fx2;
            sum2S +=fx2*fx2;
          }
          }
        double mean=(sum+sum2)/N;
        double sigma= fabs(sum-sum2)/ N * V; //sqrt(sumS/N-mean*mean);
        gsl_vector_set(vec,0,mean*V);
        gsl_vector_set(vec,1,sigma);
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
  printf("Part B\n");
  int dim3 = 3;
  double f3(int dim3, double *x3){
    return 1/((1 - cos(x3[0]) * cos(x3[1]) * cos(x3[2])) * M_PI * M_PI * M_PI);
  }
  double a3[] = {0,0,0};
  double b3[] = {M_PI,M_PI,M_PI};
  int N3 = 1e7;
  gsl_vector* vec3 = gsl_vector_alloc(dim3);

  HC_MC(dim3,f3,a3,b3,N3,vec3);
  vector_print(stdout,vec3);
  printf("∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]⁻¹ = 1.3932039296856768591842462603255\n");
  printf("Halton-Corput Monte Carlo integration with 1e7 points retuns %.10g\n",gsl_vector_get(vec3,0));
  printf("with error = %.10g\n\n", gsl_vector_get(vec3,1));

  return 0;
}
