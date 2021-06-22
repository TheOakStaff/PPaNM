#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double vector_len(gsl_vector* vec, int n);

void QR_decomp(gsl_matrix* A, gsl_matrix* R, int n, int m){
  gsl_vector *a[m], *e[m], *u[m], *b;
  for (size_t i = 0; i < m; i++) {
    a[i] = gsl_vector_alloc(n);
    gsl_matrix_get_col(a[i], A, i);
    e[i] = gsl_vector_alloc(n);
    u[i] = gsl_vector_alloc(n);
  b = gsl_vector_alloc(n);
  }
  gsl_vector_memcpy(u[0],a[0]);
  gsl_vector_memcpy(e[0],a[0]);
  gsl_vector_scale(e[0],1. / sqrt(vector_len(u[0],n)));
  int count;
  double z;
  for (size_t i = 1; i < m; i++) {
    gsl_vector_memcpy(u[i],a[i]);
    count = i;
    while (count) {
      gsl_vector_memcpy(b,e[count-1]);
      gsl_blas_ddot(a[i],e[count-1],&z);
      gsl_vector_scale(b,z);
      gsl_vector_sub(u[i],b);
      count--;
    }
    gsl_vector_memcpy(e[i],u[i]);
    gsl_vector_scale(e[i],1. / sqrt(vector_len(u[i],n)));
  }

  gsl_vector_free(b);
  for (size_t j = 0; j < m; j++) {
    gsl_matrix_set_col(A,j,e[j]);
    for (size_t i = 0; i <= j; i++) {
      gsl_blas_ddot(a[j],e[i],&z);
      gsl_matrix_set(R,i,j,z);
    }
  }
}
