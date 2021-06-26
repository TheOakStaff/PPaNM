#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// HOMEMADE INCLUSIONS //

#include "My_vec_calc.h"
#include "SVD.h"

// MAIN FUNCTION //

int main() {

	gsl_matrix* A = gsl_matrix_alloc(4,4); // (4,5) to test tall matrix
	gsl_matrix* V = gsl_matrix_alloc(A->size1,A->size2);
	gsl_matrix* Acopy = gsl_matrix_alloc(A->size1,A->size2);
	gsl_matrix* B = gsl_matrix_alloc(A->size2,A->size2);
	gsl_matrix* B2 = gsl_matrix_alloc(A->size2,A->size2);
	gsl_matrix* U = gsl_matrix_alloc(A->size1,A->size1);
	gsl_matrix* D = gsl_matrix_alloc(A->size2,A->size2);


	double S1[] = {4,-30,60,35};
  double S2[] = {-30,300,-675,420};
  double S3[] = {60,-675,1620,-1050};
  double S4[] = {-35,420,-1050,700};
//	double S5[] = {-35,420,-1050,700}; Enable to test tall matrix

  for (size_t j = 0; j < A->size2; j++) {
    gsl_matrix_set(A,0,j,S1[j]);
    gsl_matrix_set(A,1,j,S2[j]);
    gsl_matrix_set(A,2,j,S3[j]);
    gsl_matrix_set(A,3,j,S4[j]);
//		gsl_matrix_set(A,4,j,S5[j]); Enable to test tall matrix
  }

	gsl_matrix_memcpy(Acopy,A);
	printf("Starting matrix A =\n");
	matrix_print(stdout,Acopy);

	JSVD(A, V, U, D);

	printf("Final matrix A =\n");
	matrix_print(stdout,A);
	printf("Matrix V =\n");
	matrix_print(stdout,V);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,V,V,0,B);
  printf("I from V^TV\n");
  matrix_print(stdout,B);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,U,U,0,B);
  printf("I from U^TU\n");
  matrix_print(stdout,B);
	printf("D =\n");
	matrix_print(stdout,D);
	gsl_matrix_set_identity(B);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,U,D,0,B);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,B,V,0,B2);
	printf("A from UDV^T\n");
	matrix_print(stdout,B2);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(Acopy);
	gsl_matrix_free(B);
	gsl_matrix_free(B2);
	gsl_matrix_free(U);
	gsl_matrix_free(D);
	return 0;
}
