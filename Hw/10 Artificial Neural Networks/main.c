#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include "minimization.h"
#include <assert.h>


double activation(double x) {
  return x*exp(-1.0*(x*x));
}

void neuron(
  double activation(double x),
  gsl_vector *in,
  gsl_vector *out,
  gsl_vector *p,
  int neurons,
  int neuron_id)
  {
  int len = p->size;
  double x = 0;
  double r;
  for (int i = 0; i < in->size; i++) {
    x += gsl_vector_get(in, i);
  }

  r = activation((x - gsl_vector_get(p, 0 + len/neurons*neuron_id)) / gsl_vector_get(p, 1 + len/neurons*neuron_id));

  for (int i = 0; i < len/neurons - 2; i++) {  // loop over size of p's w vector
    gsl_vector_set(out, i, gsl_vector_get(p, i + len/neurons*neuron_id) * r);
  }
}

void network(
  double activation(double x),
  gsl_vector *in,
  gsl_vector *out,
  gsl_vector *p,
  int neurons)
  {
  int n = neurons;
  gsl_vector *buf[n];

  for (int i = 0; i < n; i++) {
    buf[i] = gsl_vector_alloc(out->size);
  }

  for (int i = 0; i < n; i++) {
    neuron(activation, in, buf[i], p, neurons, i);
    //vector_print(stdout, buf[i]);
  }
  double y;
  for (int i = 0; i < out->size; i++) {
    y = 0;
    for (int j = 0; j < n; j++) {
      y += gsl_vector_get(buf[j], i);
    }
    gsl_vector_set(out, i, y);
  }

  for (int i = 0; i < n; i++) {
    gsl_vector_free(buf[i]);
  }
}


double cost(            //CURRENT VERSION
  double activation(double x),
  gsl_matrix *data_in,          //maybe rewrite to take gsl_vector, ala x in code, and do organization in learn function. in order to fit with minimization algorithm
  gsl_matrix *data_out,
  gsl_vector *params,
  int neurons)
//  int neuron_id           // should be parsed to the num_grad function to select the neuron which parameters are optimized
  {
  int n = data_in->size1;
  double sum = 0;
  double sqdiv;

  gsl_vector *x = gsl_vector_alloc(data_in->size2);
  gsl_vector *y = gsl_vector_alloc(data_out->size2);
  gsl_vector *fx = gsl_vector_alloc(data_out->size2);

  for (int i = 0; i < n; i++) {
    gsl_matrix_get_row(x, data_in, i);
    gsl_matrix_get_row(y, data_out, i);

    network(activation, x, fx, params, neurons);
//    vector_print(stdout, fx);

    for (int i = 0; i < fx->size; i++) {
      sqdiv = gsl_vector_get(fx, i) - gsl_vector_get(y, i);
      sum += sqdiv*sqdiv;
    }
  }

  sum /= n;

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(fx);
//  printf("%g\n", sum);
  return sum;
}
gsl_vector* learn(
  double activation(double x),
  double cost(double activation(double x), gsl_matrix *data_in, gsl_matrix *data_out, gsl_vector *params, int neurons),
  gsl_matrix *data_in,
  gsl_matrix *data_out,
  int neurons,
  double eps)
  {

  int p_len = 2 + data_out->size2;
  gsl_vector *p = gsl_vector_alloc(p_len * neurons);
  gsl_vector_set_all(p, 3.);

  //        COST SECTION
  vector_print(stdout, p);
  qNewton(activation, cost, data_in, data_out, p, neurons, eps);
  vector_print(stdout, p);
  return p;
}

int main() {
  int neurons = 8;
  int n = 20, m_in = 1, m_out = 1;
  int seed = 5;
  double eps = 1e-3;

  gsl_matrix *data_in = gsl_matrix_alloc(n, m_in);
  gsl_matrix *data_out = gsl_matrix_alloc(n, m_out);

  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m_in; y++) {
      gsl_matrix_set(data_in, x, y, x);
      gsl_matrix_set(data_out, x, y, 2*x);
    }
  }

  gsl_vector *p;
  //matrix_print(stdout, data_in);
  //matrix_print(stdout, data_out);

  p = learn(activation, cost, data_in, data_out, neurons, eps);

  gsl_vector *x = gsl_vector_alloc(m_in);
  gsl_vector *y = gsl_vector_alloc(m_out);

  gsl_vector_set(x, 0, 5.1);

  network(activation, x, y, p, neurons);

  FILE* stream = fopen("data.txt","w");
  for (int i = 0; i < data_in->size1; i++) {
    fprintf(stream, "%g %g\n", gsl_matrix_get(data_in,i,0),gsl_matrix_get(data_out,i,0));
  }
  fclose(stream);

  FILE* stream2 = fopen("Tabdata.txt","w");
  for (double i = -2; i < 2; i+=1./8.) {
    fprintf(stream2, "%g %g\n",i,activation(i));
  }
  fclose(stream2);

  FILE* stream3 = fopen("FindPoint.txt","w");
    fprintf(stream3, "%g %g\n",gsl_vector_get(x,0),gsl_vector_get(y,0));
  fclose(stream3);


  gsl_matrix_free(data_in);
  gsl_matrix_free(data_out);
  gsl_vector_free(p);

  return 0;
}
