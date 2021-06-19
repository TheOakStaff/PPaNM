double vector_len(gsl_vector* vec){
  int n = vec->size;
  double val = 0.;
  for (size_t i = 0; i < n; i++) {
    val += gsl_vector_get(vec,i) * gsl_vector_get(vec,i);
  }
  return sqrt(val);
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R){
  int n = A->size1;
  int m = A->size2;
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
  gsl_vector_scale(e[0],1. / sqrt(vector_len(u[0])));
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
    gsl_vector_scale(e[i],1. / sqrt(vector_len(u[i])));
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

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x) {
  int n = R->size1;
  gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
  for(int i=n-1; i >= 0; i--){
    double s=gsl_vector_get(x, i);
    for(int k = i+1; k < n; k++){
      s = s - gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}

void GS_inverse(gsl_matrix *A, gsl_matrix *B){
  gsl_matrix *R = gsl_matrix_alloc(A->size2,A->size2);
  gsl_matrix *Acopy = gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(Acopy,A);
  GS_decomp(Acopy, R);
  int n = B->size1;
  gsl_matrix_set_identity(B);
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buff = gsl_vector_alloc(n);
  for(int i = 0; i < n; i++){
    gsl_matrix_get_col(e, B, i);
    GS_solve(Acopy,R,e,buff);
    gsl_matrix_set_col(B,i,buff);
  }
  gsl_vector_free(e);
  gsl_vector_free(buff);
  gsl_matrix_free(R);
}



int matrix_print(FILE* stream, const gsl_matrix *X){
  int n = X->size1;
  int m = X->size2;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.3f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n\n");
}



double vector_print(FILE* stream, gsl_vector* vec){
  int n = vec->size;
  for (size_t i = 0; i < n; i++) {
    fprintf(stream,"%g ",gsl_vector_get(vec,i));
  }
  fprintf(stream,"\n\n");
}
