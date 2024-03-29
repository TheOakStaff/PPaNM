#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <time.h>

// HOMEMADE INCLUSIONS //

#include "My_vec_calc.h"
#include "SVD.h"

// MAIN FUNCTION //

int main() {
	int n = 4;

	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* U = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);

	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* B2 = gsl_matrix_alloc(n,n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(A,i,j,(double)rand()/(double)RAND_MAX*100);
		}
	}


/* // --- Known Test matrix ---
// Making the matrix, using example matrix from:
// https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
	double S1[] = {4,-30,60,35};
  double S2[] = {-30,300,-675,420};
  double S3[] = {60,-675,1620,-1050};
  double S4[] = {-35,420,-1050,700};

  for (size_t j = 0; j < n; j++) {
    gsl_matrix_set(A,0,j,S1[j]);
    gsl_matrix_set(A,1,j,S2[j]);
    gsl_matrix_set(A,2,j,S3[j]);
    gsl_matrix_set(A,3,j,S4[j]);
  }
*/
	// RUNNING THE SVD algorithm and printing results
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
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,V,V,0,B);
	printf("I from VV^T\n");
	matrix_print(stdout,B);

	printf("Matrix U =\n");
	matrix_print(stdout,U);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,U,U,0,B);
  printf("I from U^TU\n");
	matrix_print(stdout,B);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,U,U,0,B);
	printf("I from UU^T\n");
  matrix_print(stdout,B);
	printf("D =\n");
	matrix_print(stdout,D);
	gsl_matrix_set_identity(B);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,U,D,0,B);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,B,V,0,B2);
	printf("A from UDV^T\n");
	matrix_print(stdout,B2);

// ----------------TESTING LARGE MATRIX SPEED---------------- //
	FILE* stream2 = fopen("TimeTest.txt","w");
	n = 400;
	clock_t start;
	clock_t end;
	double time_used;

	// Createing matrix and vector needed
	gsl_matrix *M = gsl_matrix_alloc(n,n);
	gsl_matrix *MV = gsl_matrix_alloc(n,n);
	gsl_matrix *MU = gsl_matrix_alloc(n,n);
	gsl_matrix *MD = gsl_matrix_alloc(n,n);

	gsl_matrix *N = gsl_matrix_alloc(n,n);
	gsl_matrix *NV = gsl_matrix_alloc(n,n);
	gsl_vector *NS = gsl_vector_alloc(n);

	gsl_matrix *BUF = gsl_matrix_alloc(n,n);

	// pseudo-Random matrix with numbers n = 0... 100)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(M,i,j,(double)rand()/(double)RAND_MAX*100);
		}
	}
	gsl_matrix_memcpy(N,M); // ensure equal starting params

	// TESTING My_SVD
	start = clock();
	JSVD(M, MV, MU, MD);
	end = clock();
	time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
	fprintf(stream2, "My_SVD function took %g seconds to compute SVD of a %dx%d matrix\n",time_used,n,n);
	// TESTING GSL one-sided Jacobi SVD
	start = clock();
	gsl_linalg_SV_decomp_jacobi(N, NV, NS);
	end = clock();
	time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
	fprintf(stream2, "GSL SVD function took %g seconds to compute SVD of a %dx%d matrix\n",time_used,n,n);

	fprintf(stream2, "\n\n\nComparison between the two resulting V^TV matrices:\n");
	fprintf(stream2,"My_SVD Matrix V^TV =\n");
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,MV,MV,0,BUF);
	matrix_print(stream2,BUF);
	fprintf(stream2,"\n");
	fprintf(stream2,"GSL SVD Matrix V^TV=\n");
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,NV,NV,0,BUF);
	matrix_print(stream2,BUF);
	fclose(stream2);

	// FREE MEMORY

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(Acopy);
	gsl_matrix_free(B);
	gsl_matrix_free(B2);
	gsl_matrix_free(U);
	gsl_matrix_free(D);
	gsl_matrix_free(M);
	gsl_matrix_free(MV);
	gsl_matrix_free(MU);
	gsl_matrix_free(MD);
	gsl_matrix_free(N);
	gsl_matrix_free(NV);
	gsl_vector_free(NS);
	gsl_matrix_free(BUF);

	return 0;
}
