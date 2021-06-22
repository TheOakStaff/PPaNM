#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"

void AtimesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta);
  double s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void jacobi_SVD(gsl_matrix *A, gsl_matrix *V){
	int n = A->size1;
	int nStep = 0;
	int breaker = 0;
	double a00;
	double a10;
	double a11;
	gsl_vector *A0 = gsl_vector_alloc(n);
	gsl_vector *A1 = gsl_vector_alloc(n);
	while (breaker == 0 && nStep < 1e7){
		nStep+=1;
		for (int i = 0; i < n-1; i++) {
			gsl_matrix_get_row(A0,A,i);
			gsl_matrix_get_row(A1,A,i+1);
			gsl_blas_ddot(A0,A0,&a00);
			gsl_blas_ddot(A1,A0,&a10);
			gsl_blas_ddot(A1,A1,&a11);
			double theta = 0.5*atan2(2*a10,a11-a00);
			double c=cos(theta),s=sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;

			AtimesJ(A,i,i+1,theta);
			if (/* condition */) {
				/* code */
			}
		}
	}
	printf("nStep =%d\n",nStep);
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
	printf("A=\n");
  matrix_print(stdout,A);
	printf("V=\n");
	matrix_print(stdout,A);
  jacobi_SVD(A,V);
	printf("A after SVD=\n");
  matrix_print(stdout,A);
	printf("V after SVD=\n");
	matrix_print(stdout,V);
/*
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,B);
	printf("A^TA=\n");
	matrix_print(stdout,B);
	*/
  return 0;
}
