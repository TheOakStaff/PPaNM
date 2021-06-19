#include "stdio.h"
#include "math.h"
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
  double delta_x = DBL_EPSILON;
  gsl_matrix *J = gsl_matrix_alloc(dim,dim);
  gsl_matrix *JR = gsl_matrix_alloc(dim,dim);
  gsl_vector *fx = gsl_vector_alloc(dim);
  gsl_vector *df = gsl_vector_alloc(dim);
  gsl_vector *Dx = gsl_vector_alloc(dim);
  gsl_vector *y = gsl_vector_alloc(dim);
  gsl_vector *fy = gsl_vector_alloc(dim);
  while(TRUE) {
    f(x,df);
    for (int j = 0; j < dim; j++) {
      gsl_vector_set(x,j,gsl_vector_get(x,j)+delta_x);
      gsl_vector_memcpy(fx,df);
      f(x,df);
      gsl_vector_sub(df,fx);
      gsl_matrix_set_col(J,j,gsl_vector_scale(df,1/delta_x));
      gsl_vector_set(x,j,gsl_vector_get(x,j) - delta_x);
      }
      GS_decomp(J,JR);
      gsl_vector_scale(fx,-1);
      GS_solve(J,JR,fx,Dx);
      gsl_vector_scale(fx,-1)
      double s = 2;
      while (TRUE) {
        s/=2.;
        gsl_vector_memcpy(y,Dx);
        gsl_vector_scale(y,s);
        gsl_vector_add(y,x);
        f(y,fy);
        if (gsl_vector_len(fy)<(1 - s / 2.) * gsl_vector_len(fx) || s<0.02) {
          break;
        }
      }
      gsl_vector_memcpy(x,y);
      gsl_vector_memcpy(fx,fy);
      if (gsl_vector_len(Dx) < delta_x || gsl_vector_len(fx) < aps) {
        break;
      }
    }
  }
}

void fun(gsl_vector *x, gsl_vector *fx){
  double GradX = -2 * (1-gsl_vector_get(x,0)) + (-2 * gsl_vector_get(x,0)) * 2 * 100 * (gsl_vector_get(x,1) - gsl_vector_get(x,0) * gsl_vector_get(x,0));
  dobule GradY = 2 * 100 * (gsl_vector_get(x,1) - gsl_vector_get(x,0) * gsl_vector_get(x,0));
  gsl_vector_set(fx,0,GradX);
  gsl_vector_set(fx,1,GradY);
}

int main() {
  /* code */
  return 0;
}
