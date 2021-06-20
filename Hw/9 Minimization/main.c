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
  double a = 2, b = 4;
  return 1 + (x-a)*(x-a) + (y-b)*(y-b);
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
      step = fabs(minimum_i) * delta_x
    }
    gsl_vector_set(minimum,i,minimum_i + step);
    gsl_vector_set(grad,i,(f(minimum)-fVal) / step);
    gsl_vector_set(minimum,i,minimum_i - step);
  }
}

void qNM(double f(gsl_vector*), gsl_vector *minimum, double esp){
  double delta_x = 2e-10;
  int dim = minimum->size;

  int nSteps = 0;
  int nReset = 0;
  int nScale = 0;

  gsl_matrix *inv_H = gsl_matrix_alloc(dim,dim), *I = gsl_matrix_alloc(dim,dim);
  gsl_matrix_set_identity = (inv_H);
  gsl_matrix_set_identity = (I);

}
