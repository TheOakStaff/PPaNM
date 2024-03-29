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

void JtimesA(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta);
  double s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}

void Jdiag(gsl_matrix* A, gsl_matrix* V){
  gsl_matrix_set_identity(V);
  int changed;
  do{
  	changed=0;
  	for(int p=0;p<A->size1-1;p++)
  	for(int q=p+1;q<A->size1;q++){
  		double apq=gsl_matrix_get(A,p,q);
  		double app=gsl_matrix_get(A,p,p);
  		double aqq=gsl_matrix_get(A,q,q);
  		double theta=0.5*atan2(2*apq,aqq-app);
  		double c=cos(theta),s=sin(theta);
  		double new_app=c*c*app-2*s*c*apq+s*s*aqq;
  		double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
  		if(new_app!=app || new_aqq!=aqq) // do rotation
  			{
  			changed=1;
  			AtimesJ(A,p,q, theta);
  			JtimesA(A,p,q,-theta); // A←J^T*A*J
  			AtimesJ(V,p,q, theta); // V←V*J
  			}
  	}
  }while(changed!=0);
}



int main() {
  int n;
  n = 4;

  gsl_matrix* A = gsl_matrix_alloc(n,n);
  gsl_matrix* V = gsl_matrix_alloc(n,n);
  gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
  gsl_matrix* Areform = gsl_matrix_alloc(n,n);
  gsl_matrix* buff = gsl_matrix_alloc(n,n);
  gsl_matrix* D = gsl_matrix_alloc(n,n);
  gsl_matrix* B = gsl_matrix_alloc(n,n);

  // https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
  double S1[] = {4,-30,60,-35};
  double S2[] = {-30,300,-675,420};
  double S3[] = {60,-675,1620,-1050};
  double S4[] = {-35,420,-1050,700};

  for (size_t j = 0; j < n; j++) {
    gsl_matrix_set(A,0,j,S1[j]);
    gsl_matrix_set(A,1,j,S2[j]);
    gsl_matrix_set(A,2,j,S3[j]);
    gsl_matrix_set(A,3,j,S4[j]);
  }

  gsl_matrix_memcpy(Acopy,A);

  printf("Starting matrix A =\n");
  matrix_print(stdout,Acopy);

  Jdiag(A, V);

  printf("Final matrix A =\n");
  matrix_print(stdout,A);

  printf("Eigen vector matrix V =\n");
  matrix_print(stdout,V);



  gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,V,Acopy,0,buff);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,buff,V,0,D);

  printf("D from V^TAV\n");
  matrix_print(stdout,D);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,V,D,0,buff);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,buff,V,0,Areform);

  printf("A from VDV^T\n");
  matrix_print(stdout,Areform);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,V,V,0,B);

  printf("I from V^TV\n");
  matrix_print(stdout,B);

  // -----------------Part B-----------------

  n=60;
  double s=1.0/(n+1);
  gsl_matrix* H = gsl_matrix_alloc(n,n);
  for(int i=0;i<n-1;i++){
    gsl_matrix_set(H,i,i,-2);
    gsl_matrix_set(H,i,i+1,1);
    gsl_matrix_set(H,i+1,i,1);
    }
  gsl_matrix_set(H,n-1,n-1,-2);
  gsl_matrix_scale(H,-1 / s / s);

  gsl_matrix *W = gsl_matrix_alloc(n,n);
  Jdiag(H,W);

  printf("k\tcalc\texact\n");
  for (int k=0; k < n/3; k++){
      double exact = M_PI*M_PI*(k+1)*(k+1);
      double calculated = gsl_matrix_get(H,k,k);
      printf("%i\t%g\t%g\n",k,calculated,exact);
  }

  printf("The calculated and exact values differ slightly and the difference becomes more pronounced with higher k values\n\n");

  double wavefunc(int j, int i, int n){
    double x;
    x = sin(((double)j+1.0) * M_PI * ((double)((i + 1.0) / (n + 1))));
    return x;
  }


  FILE* stream2 = fopen("data2.txt","w");
    fprintf(stream2,"0 0 0 0 0 0 0\n");
    for(int i=0;i<n;i++)
  	fprintf(stream2,"%g %g %g %g %g %g %g\n",(i+1.0)/(n+1), 5.52*gsl_matrix_get(W,i,0), 5.52*gsl_matrix_get(W,i,1), 5.52*gsl_matrix_get(W,i,2),wavefunc(0,i,n),-wavefunc(1,i,n),wavefunc(2,i,n));
    fprintf(stream2,"1 0 0 0 0 0 0\n");
  fclose(stream2);

  printf("The numerical and analytical results are equal if the Eigenvectors are normalized\n\n");

  return 0;
}
