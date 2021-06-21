#ifndef MY_MINIMIZER_H
#define MY_MINIMIZER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void num_grad(double f(gsl_vector*), gsl_vector* minimum, gsl_vector *grad);

void qNM(double f(gsl_vector*), gsl_vector *minimum, double esp);

#endif
