#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double vector_print(FILE* stream, gsl_vector* vec, int n){
  for (size_t i = 0; i < n; i++) {
    fprintf(stream,"%g ",gsl_vector_get(vec,i));
  }
  fprintf(stream,"\n");
}
