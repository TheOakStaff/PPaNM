#ifndef VEC_CALC_H
#define VEC_CALC_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double gsl_vector_length(gsl_vector* vec);

void matrix_print(FILE* stream, const gsl_matrix *X);

void vector_printEX(FILE* stream, gsl_vector* vec, gsl_vector *m);

void vector_print(FILE* stream, gsl_vector* vec);

void GS_decomp(gsl_matrix *A, gsl_matrix *B);

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);

void QR_inverse(gsl_matrix *Q, gsl_matrix *B);

#endif
