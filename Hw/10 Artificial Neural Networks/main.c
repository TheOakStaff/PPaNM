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

typedef struct Params{
  double a;
  double b;
  gsl_vector *w;
};

double Gauss(double x){
  return exp(-1.0*x*x);
}

void neuron(double f(double),gsl_vector *x, gsl_vector *fx, Params  *p){
  double val = 0;
  double r;

  for (int i = 0; i < x->size; i++) {
    val+= gsl_vector_get(x,i);
  }

  r = f((x - p->a) / p->b);

  for (int i = 0; i < fx->size; i++) {
    gsl_vector_set(fx,i,gsl_vector_get(p->w, i) * r);
  }
}
