#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"

double fun_k(int i, double x){
  switch (i) {
    case 0: return 1.0; break;
    case 1: return -x; break;
    default: return 0; break;
  }
}

void exp_fit(gsl_vector *x, gsl_vector *y);

int main() {
  int n = 9;
  int c_size = 2;
  gsl_vector *t = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *dy = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(c_size);
  gsl_vector *dc = gsl_vector_alloc(c_size);

  double Time[] = {1,2,3,4,6,9,10,13,15};
  double Active[] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

  for (size_t i = 0; i < n; i++) {
    gsl_vector_set(t,i,Time[i]);
    gsl_vector_set(y,i,Active[i]);
    gsl_vector_set(dy,i,Active[i]/20.);
  }

  gsl_matrix* A = gsl_matrix_alloc(t->size,c->size);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < c_size; j++) {
      gsl_matrix_set(A,i,j,fun_k(j,gsl_vector_get(t,i)/gsl_vector_get(dy,i)));
    }
  }

  printf("time(days)\n");
  vector_print(stdout,t);
  printf("Activity\n");
  vector_print(stdout,y);
  printf("Error\n");
  vector_print(stdout,dy);
  printf("Matrix A\n");
  matrix_print(stdout,A);

  return 0;
}
