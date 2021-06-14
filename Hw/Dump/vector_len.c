#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double vector_len(gsl_vector* vec, int n){
  double val = 0.;
  for (size_t i = 0; i < n; i++) {
    val += gsl_vector_get(vec,i) * gsl_vector_get(vec,i);
  }
  return val;
}
