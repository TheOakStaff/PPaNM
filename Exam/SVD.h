#ifndef MY_SVD_H
#define MY_SVD_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void AtimesJ(gsl_matrix* A, int p, int q, double theta);

void JSVD(gsl_matrix *A, gsl_matrix *V, gsl_matrix *U, gsl_matrix *D);

#endif
