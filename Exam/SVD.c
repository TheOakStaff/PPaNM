#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"
#include "SVD.h"

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

void JSVD(gsl_matrix *A, gsl_matrix *V, gsl_matrix *U, gsl_matrix *D){
  gsl_matrix_set_identity(V);
  int nSteps = 0;
  int doCheck = 1;
  gsl_vector *ap = gsl_vector_alloc(A->size1);
  gsl_vector *aq = gsl_vector_alloc(A->size1);
  gsl_vector *buf = gsl_vector_alloc(A->size1);
  gsl_vector *vec = gsl_vector_alloc(A->size1);
  double apq, app, aqq;

  while(doCheck == 1) {
    nSteps++;
    doCheck = 0;
    for (int p = 0; p < A->size2; p++) {
      for (int q = p+1; q < A->size2; q++) {

        gsl_matrix_get_col(ap, A, p);
        gsl_matrix_get_col(aq, A, q);
    		gsl_blas_ddot(ap,aq,&apq);
        gsl_blas_ddot(ap,ap,&app);
        gsl_blas_ddot(aq,aq,&aqq);

    		double theta = 0.5*atan2(2*apq,aqq-app);
    		double c=cos(theta),s=sin(theta);
    		double new_app=c*app - s*apq;
    		double new_aqq=s*apq + c*aqq;

    		if(new_app != app || new_aqq != aqq) {
    			doCheck=1;
    			AtimesJ(A,p,q, theta);
    			AtimesJ(V,p,q, theta);
    		}
    	}
    }
  }

  for (int i = 0; i < A->size1; i++) {
    gsl_matrix_get_col(vec, A, i);
    gsl_matrix_set(D,i,i,vector_len(vec));

    gsl_vector_memcpy(buf,vec);
    gsl_vector_scale(buf,1/vector_len(vec));
    gsl_matrix_set_col(U,i,buf);

//      gsl_matrix_set(U,j,i,gsl_matrix_get(A,j,i)/vector_len(vec));
  }



  gsl_vector_free(ap);
  gsl_vector_free(aq);
  gsl_vector_free(buf);
  gsl_vector_free(vec);
}
