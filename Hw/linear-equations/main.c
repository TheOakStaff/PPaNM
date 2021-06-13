#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double gsl_vector_len(gsl_vector* vec, int n){
  double val = 0.;
  for (size_t i = 0; i < n; i++) {
    val += gsl_vector_get(vec,i) * gsl_vector_get(vec,i);
  }
  printf("vector sum = %f\n",val);
  return val;
}

double gsl_vector_print(FILE* stream, gsl_vector* vec, int n){
  for (size_t i = 0; i < n; i++) {
    fprintf(stream,"%g ",gsl_vector_get(vec,i));
  }
  fprintf(stream,"\n");
}

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
  gsl_vector_scale(e[0],1. / sqrt(gsl_vector_len(u[0],n)));
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
    gsl_vector_scale(e[i],1. / sqrt(gsl_vector_len(u[i],n)));
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
  int n = 6;
  int m = 4;
  int size_m = n*m;
  int seed = 2;
  gsl_matrix *A = gsl_matrix_alloc(n,m);
  gsl_matrix *R = gsl_matrix_alloc(m,m);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(A,i,j,(double)rand_r(&seed) / (double)RAND_MAX*10);
    }
  }

  FILE* stream = fopen("data.txt","w");
  print_matrix(stream,n,m,A);
  fclose(stream);
  stream = fopen("data.txt","a");
  print_matrix(stream,m,m,R);
  fclose(stream);

  GS_decomp(A, R, n, m);

  stream = fopen("data.txt","a");
  fprintf(stream, "Q\n");
  print_matrix(stream,n,m,A);
  fclose(stream);
  stream = fopen("data.txt","a");
  fprintf(stream, "R\n");
  print_matrix(stream,m,m,R);
  fclose(stream);

  return 0;
}
