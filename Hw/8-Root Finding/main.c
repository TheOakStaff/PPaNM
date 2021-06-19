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



void newton(
  void f(gsl_vector *x, gsl_vector *fx),
  gsl_vector *x,
  double eps)
  {
  int dim = x->size;
  double delta_x = sqrt(DBL_EPSILON);
  gsl_matrix *J = gsl_matrix_alloc(dim,dim);
  gsl_matrix *JR = gsl_matrix_alloc(dim,dim);
  gsl_vector *fx = gsl_vector_alloc(dim);
  gsl_vector *df = gsl_vector_alloc(dim);
  gsl_vector *Dx = gsl_vector_alloc(dim);
  gsl_vector *y = gsl_vector_alloc(dim);
  gsl_vector *fy = gsl_vector_alloc(dim);
  gsl_vector *buf = gsl_vector_alloc(dim);
  int count = 0;
  while(1) {
    count++;
    f(x,df);
    for (int j = 0; j < dim; j++) {
      gsl_vector_memcpy(buf,df);
      gsl_vector_set(x,j,gsl_vector_get(x,j) + delta_x);
      gsl_vector_memcpy(fx,buf);
      f(x,buf);
      gsl_vector_sub(buf,fx);
      gsl_vector_scale(buf,1./delta_x);
      gsl_matrix_set_col(J,j,buf);
      gsl_vector_set(x,j,gsl_vector_get(x,j) - delta_x);
      }
    GS_decomp(J,JR);
    gsl_vector_scale(fx,-1.);
    GS_solve(J,JR,fx,Dx);
    gsl_vector_scale(fx,-1.);
    double s = 2.;
    while (1) {
      s/=2.;
      gsl_vector_memcpy(y,Dx);
      gsl_vector_scale(y,s);
      gsl_vector_add(y,x);
      f(y,fy);
      if (vector_len(fy)<(1 - s / 2.) * vector_len(fx) || s<0.02) {
        break;
      }
    }
    gsl_vector_memcpy(x,y);
    gsl_vector_memcpy(fx,fy);
    if (vector_len(Dx) < delta_x || vector_len(fx) < eps  || count > 1e8) {
      break;
    }
  }
  gsl_vector_free(fx);
  gsl_vector_free(df);
  gsl_vector_free(Dx);
  gsl_vector_free(y);
  gsl_vector_free(fy);
  gsl_vector_free(buf);
  gsl_matrix_free(J);
  gsl_matrix_free(JR);
}


void fun(gsl_vector *x, gsl_vector *fx){
  double a = 1;
  double b = 100;
  double x1 = gsl_vector_get(x,0);
  double y1 = gsl_vector_get(x,1);
  double GradX = -2*a + 4*b*x1*x1*x1 - 4*b*x1*y1 + 2*x1;
  double GradY =  2*b*(y1-x1*x1);
//  double GradX = -2 * (1-gsl_vector_get(x,0)) + (-2 * gsl_vector_get(x,0)) * 2 * 100 * (gsl_vector_get(x,1) - gsl_vector_get(x,0) * gsl_vector_get(x,0));
//  double GradY = 2 * 100 * (gsl_vector_get(x,1) - gsl_vector_get(x,0) * gsl_vector_get(x,0));
  gsl_vector_set(fx,0,GradX);
  gsl_vector_set(fx,1,GradY);
}


void fun1D(gsl_vector *x, gsl_vector *fx){
  gsl_vector_set(fx,0,gsl_vector_get(x,0)*2+1);
}

void fun2D(gsl_vector *x, gsl_vector *fx){
  gsl_vector_set(fx,0,gsl_vector_get(x,0)*2);
  gsl_vector_set(fx,1,gsl_vector_get(x,1)*2);
}


int main() {

  int n = 1;
  double tolerance = 1e-8;
  gsl_vector *X_val = gsl_vector_alloc(n);
  gsl_vector_set(X_val,0,2);
  newton(fun1D,X_val,tolerance);
  vector_print(stdout,X_val);
  gsl_vector_free(X_val);


  n = 2;
  gsl_vector *X_val2 = gsl_vector_alloc(n);
  gsl_vector_set(X_val2,0,2);
  gsl_vector_set(X_val2,1,2);
  newton(fun2D,X_val2,tolerance);
  vector_print(stdout,X_val2);
  gsl_vector_free(X_val2);

  n = 2;
  gsl_vector *X_val3 = gsl_vector_alloc(n);
  gsl_vector_set(X_val3,0,0);
  gsl_vector_set(X_val3,1,3);
  printf("The initial value is (x, y) = (%g, %g)\n",gsl_vector_get(X_val3,0), gsl_vector_get(X_val3,1));
  printf("The minimum should be at (x, y)= (1, 1)\n");
  newton(fun,X_val3,tolerance);
  printf("The minimum is found at (x, y) = (%g, %g),\n\n", gsl_vector_get(X_val3,0), gsl_vector_get(X_val3,1));
  gsl_vector_free(X_val3);
  printf("The final value differs from the actual value. However, plotting the");
  printf("function using Wolfram-Alpha[1] reveals that the 'valley floor' is rather extencive.\n");
  printf("As such, our value of likely lies in an area where the \nrate of change is too low for our function to handle\n\n");
  printf("[1] https://bit.ly/3qaU34S\n\n");
  return 0;
}
