#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"

void AJ(gsl_matrix *A){

}

void jacobi_SVD(gsl_matrix *A, gsl_matrix *V){

	int n = A->size1;
	double a00;
	double a10;
	double a11;
	gsl_vector *A0;
	gsl_vector *A1;
	for (int i = 0; i < n-1; i++) {
		gsl_matrix_get_row(A0,A,i);
		gsl_matrix_get_row(A1,A,i+1);
		gsl_blas_ddot(A0,A0,&a00);
		gsl_blas_ddot(A1,A0,&a10);
		gsl_blas_ddot(A1,A1,&a11);
		double theta = atan2(2*a10,a11-a00);
	}



}











int main() {
  // USING WIKIPEDIA TEST MATRIX
  gsl_matrix *A = gsl_matrix_alloc(2,2);
  gsl_matrix *V = gsl_matrix_alloc(2,2);
	gsl_matrix *B = gsl_matrix_alloc(2,2);

  double S1[] = {2,4};
  double S2[] = {6,8};
//  double S3[] = {0,0,0,0};
//  double S4[] = {0,2,0,0};

  for (size_t j = 0; j < 2; j++) {
    gsl_matrix_set(A,0,j,S1[j]);
    gsl_matrix_set(A,1,j,S2[j]);
//    gsl_matrix_set(A,2,j,S3[j]);
//    gsl_matrix_set(A,3,j,S4[j]);
  }

  //jacobi_SVD(A,V);
	printf("A=\n");
  matrix_print(stdout,A);
	printf("V=\n");
	matrix_print(stdout,V);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,B);
	printf("A^TA=\n");
	matrix_print(stdout,B);
  return 0;
}
