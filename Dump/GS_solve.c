#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x, int n) {
  gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(x, i);
    for(int k = i+1; k < n; k++){
      s = s - gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}
