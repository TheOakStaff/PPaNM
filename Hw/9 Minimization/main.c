#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "float.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"
#include <assert.h>

double testFunction(gsl_vector *vec){
  double x = gsl_vector_get(vec,0), y = gsl_vector_get(vec,1);
  double a = -2, b = -2;
  return (x+a)*(x+a) + (y+a)*(y+a);
}

double RosValley(gsl_vector *vec){
  double x = gsl_vector_get(vec,0), y = gsl_vector_get(vec,1);
  double a = 1, b = 100;
  return (a-x)*(a-x) + b * (y-x*x)*(y-x*x);
}

double Himmelblau(gsl_vector *vec){
  double x = gsl_vector_get(vec,0), y = gsl_vector_get(vec,1);
  double a = 11, b = 7;
  return (x*x + y - a)*(x*x + y - a) + (x + y*y - b)*(x + y*y - b);
}

void num_grad(double f(gsl_vector*), gsl_vector* minimum, gsl_vector *grad){
  double delta_x = 2e-10;

  double fVal = f(minimum);
  int dim = minimum->size;
  double step;

  for (int i = 0; i < dim; i++) {

    double minimum_i = gsl_vector_get(minimum,i);

    if (fabs(minimum_i) < delta_x) {
      step = delta_x;
    }
    else{
      step = fabs(minimum_i) * delta_x;
    }
    gsl_vector_set(minimum,i,minimum_i + step);
    gsl_vector_set(grad,i,(f(minimum)-fVal) / step);
    gsl_vector_set(minimum,i,minimum_i - step);
  }
}

void qNM(double f(gsl_vector*), gsl_vector *minimum, double esp){
  double delta_x = 1e-10;
  int dim = minimum->size;

  int nSteps = 0;
  int nResets = 0;
  int nScales = 0;

  gsl_matrix *inv_H = gsl_matrix_alloc(dim,dim);
  gsl_matrix *I = gsl_matrix_alloc(dim,dim);
  gsl_matrix_set_identity(inv_H);
  gsl_matrix_set_identity(I);

  gsl_vector *grad = gsl_vector_alloc(dim);
  gsl_vector *newtStep = gsl_vector_alloc(dim);
  gsl_vector *minimumNext = gsl_vector_alloc(dim);
  gsl_vector *gradNext = gsl_vector_alloc(dim);
  gsl_vector *sol = gsl_vector_alloc(dim);
  gsl_vector *solDif = gsl_vector_alloc(dim);

  gsl_matrix* solDifsolDifT = gsl_matrix_calloc(dim, dim);

  num_grad(f,minimum,grad);
  double fVal = f(minimum);
  double fValNext;
  double scale;
  double sTg;
  double solDifTsol;

  while (nSteps < 1e5) {
    nSteps++;

    gsl_blas_dgemv(CblasNoTrans,-1,inv_H,grad,0,newtStep);
    if (vector_len(newtStep) < delta_x * vector_len(minimum) ) {
      fprintf(stderr, "Error 1: NewtonStep smaller than allowed minimum\n");
      break;
    }

    if (vector_len(grad) < esp) {
      fprintf(stderr, "Error 2: Gradient smaller than tolerance\n");
      break;
    }
    scale = 1.;

    while (1) {
      gsl_vector_memcpy(minimumNext,minimum);
      gsl_vector_add(minimumNext,newtStep);

      fValNext = f(minimumNext);

      gsl_blas_ddot(newtStep,grad,&sTg);

      if (fValNext < fVal + 0.01 *sTg) {
        nScales++;
        break;
      }

      if (scale < delta_x) {
        nResets++;
        gsl_matrix_set_identity(inv_H);
        break;
      }
      scale*=0.5;
      gsl_vector_scale(newtStep,0.5);
    }

    num_grad(f,minimumNext,gradNext);

    gsl_vector_memcpy(sol,gradNext);
    gsl_blas_daxpy(-1,grad,sol);
    gsl_vector_memcpy(solDif,newtStep);
    gsl_blas_dgemv(CblasNoTrans,-1,inv_H,sol,1,solDif);

    gsl_blas_dsyr(CblasUpper, 1.0, solDif, solDifsolDifT);
    gsl_blas_ddot(solDif,sol,&solDifTsol);

    if (fabs(solDifTsol) < 1e-12) {
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0/solDifTsol,solDifsolDifT,I,1.0,inv_H);
    }

    gsl_vector_memcpy(minimum,minimumNext);
    gsl_vector_memcpy(grad,gradNext);
    fVal = fValNext;
  }
  gsl_matrix_free(inv_H);
  gsl_matrix_free(I);
  gsl_matrix_free(solDifsolDifT);
  gsl_vector_free(grad);
  gsl_vector_free(newtStep);
  gsl_vector_free(minimumNext);
  gsl_vector_free(gradNext);
  gsl_vector_free(sol);
  gsl_vector_free(solDif);

  printf("\nnSteps = %i \nnScales = %i \nnResets = %i \nfVal = %g\n\n",nSteps,nScales,nResets,fVal);
}


int main() {
  int dim = 2;
  double tol = 1e-6;

  gsl_vector *minimum = gsl_vector_alloc(dim);
  gsl_vector_set(minimum,0,-20);
  gsl_vector_set(minimum,1,-20);
  printf("Testing using f(x,y) = (x-2)² + (y-2)² \n");
  printf("Initial values = (%g,%g)\nActual minimum located at (2,2)\n",gsl_vector_get(minimum,0),gsl_vector_get(minimum,1));
  qNM(testFunction,minimum,tol);
  printf("The found minimum is (%g,%g)\n\n", gsl_vector_get(minimum,0),gsl_vector_get(minimum,1));

  int dim2 = 2;
  double tol2 = 1e-6;

  gsl_vector *minimum2 = gsl_vector_alloc(dim2);
  gsl_vector_set(minimum2,0,2);
  gsl_vector_set(minimum2,1,2);
  printf("Testing Rosenbrocks valley function a=1, b=100 \n");
  printf("Initial values = (%g,%g)\nActual minimum located at (1,1)\n",gsl_vector_get(minimum2,0),gsl_vector_get(minimum2,1));
  qNM(RosValley,minimum2,tol2);
  printf("The found minimum is (%g,%g)\n\n", gsl_vector_get(minimum2,0),gsl_vector_get(minimum2,1));

  int dim3 = 2;
  double tol3 = 1e-6;

  gsl_vector *minimum3 = gsl_vector_alloc(dim3);
  gsl_vector_set(minimum3,0,4);
  gsl_vector_set(minimum3,1,4);
  printf("Testing Himmelblau function\n");
  printf("Initial values = (%g,%g)\nActual minimum located at \n(3,2)  (-2.80511,3.1313)  (-3.7793:-3.2831)  (3.5844,-1.84812)\n",gsl_vector_get(minimum3,0),gsl_vector_get(minimum3,1));
  qNM(Himmelblau,minimum3,tol3);
  printf("The found minimum is (%g,%g)\n", gsl_vector_get(minimum3,0),gsl_vector_get(minimum3,1));
  return 0;
}
