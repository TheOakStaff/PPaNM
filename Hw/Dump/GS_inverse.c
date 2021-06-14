#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x, int n);

void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B,gsl_matrix *X, int n){
  gsl_matrix_set_identity(B);
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buff = gsl_vector_alloc(n);
  for(int i = 0; i < n; i++){
    gsl_matrix_get_col(e, B, i);
    GS_solve(Q,R,e,buff,n);
    gsl_matrix_set_col(X,i,buff);
  }
  gsl_vector_free(e);
  gsl_vector_free(buff);
}
