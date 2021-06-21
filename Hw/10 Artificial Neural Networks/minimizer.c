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
