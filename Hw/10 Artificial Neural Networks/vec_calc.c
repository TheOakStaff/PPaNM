#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include <assert.h>

double gsl_vector_length(gsl_vector* vec){
  int n = vec->size;
  double sum = 0.;
  for (int x = 0; x < n; x++) {
    sum += pow(gsl_vector_get(vec, x),2);
  }
  return sqrt(sum);
}

void matrix_print(FILE* stream, const gsl_matrix *X){
  int n = X->size1;
  int m = X->size2;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

void vector_printEX(FILE* stream, gsl_vector* vec, gsl_vector *m){
  int n = vec->size;
  if (!gsl_vector_isnull(m)) {
    for (int i = 0; i < m->size; i++) {
      fprintf(stream, "%g  ",gsl_vector_get(m, i));
    }
  }
  for (int x = 0; x < n; x++) {
    fprintf(stream, "%g  ",gsl_vector_get(vec, x));
  }
  fprintf(stream,"\n");
}

void vector_print(FILE* stream, gsl_vector* vec){
  int n = vec->size;
  for (int x = 0; x < n; x++) {
    fprintf(stream, "%g ",gsl_vector_get(vec, x));
  }
  fprintf(stream,"\n\n");
}

void GS_decomp(gsl_matrix *A, gsl_matrix *B){
  int n = A->size1;
  int m = A->size2;

  gsl_vector* a[m], *e[m], *u[m];
  for (size_t x = 0; x < m; x++) {
    a[x] = gsl_vector_alloc(n);
    gsl_matrix_get_col(a[x], A, x);
    u[x] = gsl_vector_alloc(n);
    e[x] = gsl_vector_alloc(n);
  }
  gsl_vector_memcpy(u[0],a[0]);
  gsl_vector_memcpy(e[0],u[0]);
  gsl_vector_scale(e[0],1./sqrt(gsl_vector_length(u[0])));
  gsl_vector *buf = gsl_vector_alloc(n);
  int level;
  double mul;
  for (int x = 1; x < m; x++) {
    level = x;
    gsl_vector_memcpy(u[x],a[x]);
    while (level) {
      gsl_vector_memcpy(buf,e[level-1]);
      gsl_blas_ddot(a[x],e[level-1],&mul);
      gsl_vector_scale(buf,mul);
      gsl_vector_sub(u[x],buf);
      level--;
    }
    gsl_vector_memcpy(e[x],u[x]);
    gsl_vector_scale(e[x],1./sqrt(gsl_vector_length(e[x])));
  }
  gsl_vector_free(buf);
  double val;
  for (int x = 0; x < m; x++) {
    gsl_matrix_set_col(A, x, e[x]);
    for (int y = 0; y <= x; y++) {
      gsl_blas_ddot(a[x], e[y], &val);
      gsl_matrix_set(B, y, x, val);
    }
  }
  for (int x = 0; x < m; x++) {
    gsl_vector_free(a[x]);
    gsl_vector_free(u[x]);
    gsl_vector_free(e[x]);
  }
}

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x) {
  int n = R->size1;

  gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(x, i);
    for(int k = i+1; k < n; k++){
      s -= gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}

void QR_inverse(gsl_matrix *Q, gsl_matrix *B){
  int n = Q->size2;
  gsl_matrix *R = gsl_matrix_alloc(n, n);
  gsl_matrix *A = gsl_matrix_alloc(Q->size1, Q->size2);
  gsl_matrix_memcpy(A, Q);
  GS_decomp(Q, R);

  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buf = gsl_vector_alloc(n);
  gsl_matrix_set_identity(B);
  for(int i = 0; i<n ; i++){
    gsl_matrix_get_col(e, B, i);
    GS_solve(Q, R, e, buf);
    gsl_matrix_set_col(B, i, buf);
  }
  gsl_matrix_memcpy(Q, A);
  gsl_vector_free(e);
  gsl_vector_free(buf);
  gsl_matrix_free(R);
  gsl_matrix_free(A);
}
