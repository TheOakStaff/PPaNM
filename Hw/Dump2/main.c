#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "float.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"
#include <assert.h>

FILE* stream = fopen("data.txt","w");
for (int i = 0; i < n; i++) {
  fprintf(stream, "%g %g\n", gsl_matrix_get(data_in,0,i),gsl_matrix_get(data_out,0,i));
}
fclose(stream);
