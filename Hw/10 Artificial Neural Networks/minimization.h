#ifndef MY_MINIMIZER_H
#define MY_MINIMIZER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void num_grad(
  double activation(double x),
  double f(double activation(double x), gsl_matrix *data_in, gsl_matrix *data_out, gsl_vector *params, int neurons),
  gsl_matrix *data_in,
  gsl_matrix *data_out,
  gsl_vector *x,
  int neurons,
  gsl_vector *grad);

void qNewton(
  double activation(double x),
  double f(double activation(double x), gsl_matrix *data_in, gsl_matrix *data_out, gsl_vector *params, int neurons),
  gsl_matrix *data_in,
  gsl_matrix *data_out,
  gsl_vector *x,
  int neurons,
  double eps);

  #endif
