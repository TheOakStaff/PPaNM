#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R, int n, int m){
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
  gsl_vector_scale(e[0],1. / sqrt(gsl_vector_sum(e[0])));
  int count;
  double z;
  for (size_t i = 0; i < m; i++) {
    count = 0;
    while (count < i) {
      sl_vector_memcpy(b,e[count]);
      gsl_blas_dsdot(a[i]),e[count],&z);
      gsl_vector_sub(u[i],gsl_vector_scale(b,z));
      count++;
    }
    gsl_vector_memcpy(e[i],u[i]);
    gsl_vector_scale(e[i],1. / sqrt(gsl_vector_sum(u[i]));
  }
  gsl_vector_free(b);
  for (size_t j = 0; j < m; j++) {
    gsl_matrix_set_col(A,i,e[i]);
    for (size_t i = 0; i < n; i++) {
      gsl_matrix_set(R,i,j,gsl_blas_dsdot(a[i]),e[j],&z);
    }
  }
}


int print_matrix(FILE* stream, int n, int m, const gsl_matrix *X){
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

int main(int argc, char const *argv[]) {
  int n = 4;
  int m = 3;
  int size_m = n*m;
  int seed = 42;
  gsl_matrix *A = gsl_matrix_alloc(n,m);
  gsl_matrix *R = gsl_matrix_alloc(n,m);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(A,i,j,rand_r(&seed) % 10);
      gsl_matrix_set(R,i,j,0);
    }
  }
  FILE* stream = fopen("data.txt","w");
  print_matrix(stream,n,m,A);
  fclose(stream);
  stream = fopen("data.txt","a");
  print_matrix(stream,n,m,R);
  fclose(stream);

  GS_decomp(A, R, n, m);

  stream = fopen("Q.txt","w");
  print_matrix(stream,n,m,A);
  fclose(stream);
  stream = fopen("QR.txt","w");
  print_matrix(stream,n,m,R);
  fclose(stream);

  return 0;
}
