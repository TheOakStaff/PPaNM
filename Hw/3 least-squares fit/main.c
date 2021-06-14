#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"

double fun_k(int i, double x){
  switch (i) {
    case 0: return 1.0; break;
    case 1: return -x; break;
    default: return 0; break;
  }
}

void exp_fit(gsl_vector *x, gsl_vector *y, gsl_vector *dy, gsl_matrix *A, gsl_vector *c, gsl_vector *dc){
  gsl_vector *b = gsl_vector_alloc(y->size);
  for (size_t i = 0; i <A->size1; i++) {
    gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
    for (size_t j = 0; j <A->size2; j++) {
      gsl_matrix_set(A,i,j,fun_k(j,gsl_vector_get(x,i))/gsl_vector_get(dy,i));
    }
  }

  gsl_matrix *Q = gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix *R = gsl_matrix_alloc(A->size2,A->size2);
  gsl_matrix_memcpy(Q,A);
  GS_decomp(Q,R);
  GS_solve(Q,R,b,c);

  gsl_matrix *INV = gsl_matrix_alloc(R->size1,R->size2);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,R);
  GS_inverse(R,INV);
  for (size_t i = 0; i <INV->size1; i++) {
      gsl_vector_set(dc,i,sqrt(gsl_matrix_get(INV,i,i)));
  }


}

int main() {
  int n = 9;
  int c_size = 2;
  gsl_vector *t = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *dy = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(c_size);
  gsl_vector *dc = gsl_vector_alloc(c_size);
  double Time[] = {1,2,3,4,6,9,10,13,15};
  double Active[] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

  for (size_t i = 0; i < n; i++) {
    gsl_vector_set(t,i,Time[i]);
    gsl_vector_set(y,i,logf(Active[i]));
    gsl_vector_set(dy,i,Active[i]/20./Active[i]);
  }

  gsl_matrix* A = gsl_matrix_alloc(t->size,c->size);

  printf("time(days)\n");
  vector_print(stdout,t);
  printf("Activity\n");
  vector_print(stdout,y);
  printf("dy (Error)\n");
  vector_print(stdout,dy);

  exp_fit(t,y,dy,A,c,dc);

  printf("Matrix A\n");
  matrix_print(stdout,A);
  printf("vector c (a,lambda)\n");
  vector_print(stdout,c);
  printf("with uncertainties:\nda = %f\t dlambda = %f\n",gsl_vector_get(dc,0),gsl_vector_get(dc,1));

  FILE* stream = fopen("fit_data.txt","w");
  for (double i = 0; i < 20; i+=0.1) {
    fprintf(stream, "%g %g %g %g\n",i,gsl_vector_get(c,0)-gsl_vector_get(c,1)*i,
    (gsl_vector_get(c,0)-gsl_vector_get(dc,0))-(gsl_vector_get(c,1)-gsl_vector_get(dc,1))*i,
    (gsl_vector_get(c,0)+gsl_vector_get(dc,0))-(gsl_vector_get(c,1)+gsl_vector_get(dc,1))*i);
  }
  fclose(stream);

  FILE* stream2 = fopen("data.txt","w");
  for (size_t i = 0; i <y->size; i++) {
    fprintf(stream2, "%g %g %g\n",gsl_vector_get(t,i),gsl_vector_get(y,i),gsl_vector_get(dy,i));
  }
  fclose(stream2);
  double halflife;
  halflife = (gsl_vector_get(c,0)-logf(117/2.))/gsl_vector_get(c,1);
  printf("Half life can be found using log(2)/lambda\n");
  printf("Using this we find the half-life to be %f days \n",log(2)/gsl_vector_get(c,1));
  printf("within the uncertainties\nLower-half-life = %f\tUpper-half-life = %f\n",log(2)/(gsl_vector_get(c,1)+gsl_vector_get(dc,1)),log(2)/(gsl_vector_get(c,1)-gsl_vector_get(dc,1)));
  printf("Ra224 has a actual half-life of 3.6319 days, so even with uncertainties the result is still of by almost half a day\n");

  return 0;
}
