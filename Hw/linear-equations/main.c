#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double gsl_vector_len(gsl_vector* vec, int n){
  double val = 0.;
  for (size_t i = 0; i < n; i++) {
    val += gsl_vector_get(vec,i) * gsl_vector_get(vec,i);
  }
  return val;
}

double gsl_vector_print(FILE* stream, gsl_vector* vec, int n){
  for (size_t i = 0; i < n; i++) {
    fprintf(stream,"%g ",gsl_vector_get(vec,i));
  }
  fprintf(stream,"\n");
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R, int n, int m){
  gsl_vector *a[m], *e[m], *u[m], *b;
  for (size_t i = 0; i < m; i++) {
    a[i] = gsl_vector_alloc(n);
    gsl_matrix_get_col(a[i], A, i);
    e[i] = gsl_vector_alloc(n);
    u[i] = gsl_vector_alloc(n);
  b = gsl_vector_alloc(n);
  }
  gsl_vector_memcpy(u[0],a[0]);
  gsl_vector_memcpy(e[0],a[0]);
  gsl_vector_scale(e[0],1. / sqrt(gsl_vector_len(u[0],n)));
  int count;
  double z;
  for (size_t i = 1; i < m; i++) {
    gsl_vector_memcpy(u[i],a[i]);
    count = i;
    while (count) {
      gsl_vector_memcpy(b,e[count-1]);
      gsl_blas_ddot(a[i],e[count-1],&z);
      gsl_vector_scale(b,z);
      gsl_vector_sub(u[i],b);
      count--;
    }
    gsl_vector_memcpy(e[i],u[i]);
    gsl_vector_scale(e[i],1. / sqrt(gsl_vector_len(u[i],n)));
  }

  gsl_vector_free(b);
  for (size_t j = 0; j < m; j++) {
    gsl_matrix_set_col(A,j,e[j]);
    for (size_t i = 0; i <= j; i++) {
      gsl_blas_ddot(a[j],e[i],&z);
      gsl_matrix_set(R,i,j,z);
    }
  }
}


int print_matrix(FILE* stream, int n, int m, const gsl_matrix *X){
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x, int n) {
  gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(x, i);
    for(int k = i+1; k < n; k++){
      s = s - gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}

/*
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x, int l){
gsl_vector* Qb = gsl_vector_alloc(l);
gsl_blas_dgemv(CblasTrans,1,Q,b,0,Qb);
double z = gsl_vector_get(Qb,l-1);
gsl_vector_set(x,l-1,z);
double val = 0.;
for (size_t i = l-1; i > 0; i--) {
  int current = i;
  while (i<l) {
    val = gsl_vector_get(Qb,i-1) - gsl_vector_get(x,i);
    gsl_vector_set(x,l-1,val);
*/

void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B,gsl_matrix *X, int n){
  gsl_matrix_set_identity(B);
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buff = gsl_vector_alloc(n);
  for(int i = 0; i < n; i++){
    gsl_matrix_get_col(e, B, i);
    GS_solve(Q,R,e,buff,n);
    gsl_matrix_set_col(X,i,buff);
  }
  gsl_vector_free(e);
  gsl_vector_free(buff);
}

int main() {
  int n = 6;
  int m = 4;
  int size_m = n*m;
  int seed = 3;
  gsl_matrix *A = gsl_matrix_alloc(n,m);
  gsl_matrix *R = gsl_matrix_alloc(m,m);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      gsl_matrix_set(A,i,j,(double)rand_r(&seed) / (double)RAND_MAX*10);
    }
  }

  int l = 4;
  gsl_vector *b = gsl_vector_alloc(l);
  gsl_vector *x = gsl_vector_alloc(l);
  gsl_matrix *L = gsl_matrix_alloc(l,l);
  gsl_matrix *K = gsl_matrix_alloc(l,l);
  for (size_t i = 0; i < l; i++) {
    gsl_vector_set(b,i,(double)rand_r(&seed) / (double)RAND_MAX*10);
    for (size_t j = 0; j < l; j++) {
      gsl_matrix_set(L,i,j,(double)rand_r(&seed) / (double)RAND_MAX*10);
    }

  }

  FILE* stream2 = fopen("data2.txt","w");
  fprintf(stream2, "vector b\n");
  gsl_vector_print(stream2,b,l);
  fprintf(stream2, "\n");
  fclose(stream2);

  FILE* stream = fopen("data.txt","w");
  fprintf(stream, "A\n");
  print_matrix(stream,n,m,A);
  fclose(stream);

  stream = fopen("data.txt","a");
  fprintf(stream, "Blank Matrix\n");
  print_matrix(stream,m,m,R);
  fclose(stream);

  GS_decomp(A, R, n, m);
  gsl_matrix *C = gsl_matrix_alloc(m,m);
  gsl_matrix *Anew = gsl_matrix_alloc(n,m);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,C);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,R,0,Anew);

  gsl_matrix *Qnew = gsl_matrix_alloc(l,l);
  gsl_matrix_memcpy(Qnew,L);

  stream2 = fopen("data2.txt","a");
  fprintf(stream2, "A\n");
  print_matrix(stream2,l,l,L);
  fclose(stream2);

  GS_decomp(L, K, l, l);
  stream2 = fopen("data2.txt","a");
  fprintf(stream, "Q\n");
  print_matrix(stream,l,l,L);
  fprintf(stream, "R\n");
  print_matrix(stream,l,l,K);
  fclose(stream2);
  GS_solve(L,K,b,x,l);

  stream2 = fopen("data2.txt","a");
  fprintf(stream2, "After back-substitution\n");
  fprintf(stream2, "Q\n");
  print_matrix(stream2,l,l,L);
  fprintf(stream2, "R\n");
  print_matrix(stream2,l,l,K);
  fprintf(stream2, "b\n");
  gsl_vector_print(stream2,b,l);
  fprintf(stream2, "\n");
  fprintf(stream2, "x\n");
  gsl_vector_print(stream2,x,l);
  fprintf(stream2, "\n");


  gsl_vector *bnew = gsl_vector_alloc(l);
  gsl_blas_dgemv(CblasNoTrans,1,Qnew,x,0,bnew);

  fprintf(stream2, "bnew\n");
  gsl_vector_print(stream2,bnew,l);
  fprintf(stream2, "\n");
  fclose(stream2);

  stream = fopen("data.txt","a");
  fprintf(stream, "Q\n");
  print_matrix(stream,n,m,A);

  fprintf(stream, "R\n");
  print_matrix(stream,m,m,R);

  fprintf(stream, "Q^TQ\n");
  print_matrix(stream,m,m,C);

  fprintf(stream, "Anew\n");
  print_matrix(stream,n,m,Anew);
  fclose(stream);

  printf("The results for A.1 can be found in data.txt\nThe results for A.2 can be found in data2.txt\n");

  // PART B

  int h = 4;
  gsl_matrix *Cnew = gsl_matrix_alloc(h,h);
  gsl_matrix *CnewCopy = gsl_matrix_alloc(h,h);
  gsl_matrix *D = gsl_matrix_alloc(h,h);
  gsl_matrix *E = gsl_matrix_alloc(h,h);
  gsl_matrix *X = gsl_matrix_alloc(h,h);
  gsl_matrix *I = gsl_matrix_alloc(h,h);
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < h; j++) {
      gsl_matrix_set(Cnew,i,j,(double)rand_r(&seed) / (double)RAND_MAX*10);
    }
  }

  gsl_matrix_memcpy(CnewCopy,Cnew);

  FILE* stream3 = fopen("data3.txt","w");

  fprintf(stream3, "Q\n");
  print_matrix(stream3,h,h,Cnew);
  fprintf(stream3, "R\n");
  print_matrix(stream3,h,h,D);

  GS_decomp(Cnew, D, h, h);
  GS_inverse(Cnew, D, E, X, h);

  fprintf(stream3, "Q after \n");
  print_matrix(stream3,h,h,Cnew);
  fprintf(stream3, "R after \n");
  print_matrix(stream3,h,h,D);

  fprintf(stream3, "B^-1\n");
  print_matrix(stream3,h,h,X);

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,CnewCopy,0,I);
  fprintf(stream3, "B^-1 * A = I\n");
  print_matrix(stream3,h,h,I);

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,CnewCopy,X,0,I);
  fprintf(stream3, "A * B^-1 = I\n");
  print_matrix(stream3,h,h,I);

  fclose(stream3);
  printf("The results for B can be found in data3.txt\n");
  return 0;
}
