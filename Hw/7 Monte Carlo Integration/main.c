#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"


void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N, gsl_vector* A){
        double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0,sum2=0,x[dim];
        for(int i=0;i<N;i++){
                for(int i=0;i<dim;i++)x[i]=a[i]+RANDOM*(b[i]-a[i]);
                double fx=f(dim,x); sum+=fx; sum2+=fx*fx;
                }
        double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
        complex result=mean*V+I*sigma*V/sqrt(N);
        return result;
}
