#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int main() {
  double a[] = {
    6.13,-2.90,5.86,
    8.08,-6.31,-3.89,
    -4.36,1.00,0.19
  };

  double x[] = {
    0,0,0
  };

  double b[] = {
    6.23, 5.37, 2.29
  };

  gsl_matrix_view A = gsl_matrix_view_array(a,3,3);

  gsl_vector_view B = gsl_vector_view_array(b,3);
  gsl_vector_view X = gsl_vector_view_array(x,3);

  gsl_linalg_HH_solve(&A.matrix, &B.vector, &X.vector);

  printf("x0= %f\nx1= %f\nx2=  %f\n",x[0],x[1],x[2]);

  double a2[] = {
    6.13,-2.90,5.86,
    8.08,-6.31,-3.89,
    -4.36,1.00,0.19
  };

  double x2[] = {
    x[0],x[1],x[2]
  };

  double b2[] = {
    0,0,0
  };

  gsl_matrix_view A2 = gsl_matrix_view_array(a2,3,3);
  gsl_vector_view X2 = gsl_vector_view_array(x2,3);
  gsl_vector_view B2 = gsl_vector_view_array(b2,3);

  gsl_blas_dgemv(CblasNoTrans,1,&A2.matrix,&X2.vector,0,&B2.vector);

  printf("\nb0=  %f\nb1=  %f\nb2=  %f\n",b2[0],b2[1],b2[2]);

  return 0;
}
