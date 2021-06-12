#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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
  return 0;
}
