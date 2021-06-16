#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "My_vec_calc.h"
#include "assert.h"

void rkstep12(
	void (*f)(int n, double t,gsl_vector *y,gsl_vector *dydt), /* the f from dy/dt=f(t,y) */
  double t,              /* the current value of the variable */
	gsl_vector *yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector *yh,             /* output: y(t+h) */
	gsl_vector *err             /* output: error estimate */
) {
  int n = yt->size;
  int i;
  gsl_vector *k0 = gsl_vector_alloc(n);
  gsl_vector *k12 = gsl_vector_alloc(n);
  gsl_vector *temp_yt = gsl_vector_alloc(n);

  f(n,t,yt,k0);
  for (i = 0; i < n; i++) {
    gsl_vector_set(temp_yt,i,gsl_vector_get(yt,i) + gsl_vector_get(k0,i) * h / 2);
  }
  f(n,t+h/2,temp_yt,k12);
  for (i = 0; i < n; i++) {
    gsl_vector_set(yh,i,gsl_vector_get(yt,i) + gsl_vector_get(k12,i) * h);
  }
  for (i = 0; i < n; i++) {
    gsl_vector_set(err,i,(gsl_vector_get(k0,i) - gsl_vector_get(k12,i)) * h / 2);
  }
}

void rkstep12SIR(
	void (*ftc)(int n, double t,gsl_vector *y,gsl_vector *dydt, double Tc, double Tr, int N), /* the f from dy/dt=f(t,y) */
  double t,              /* the current value of the variable */
	gsl_vector *yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector *yh,             /* output: y(t+h) */
	gsl_vector *err,             /* output: error estimate */
  double Tc,
  double Tr,
  int N
) {
  int n = yt->size;
  int i;
  gsl_vector *k0 = gsl_vector_alloc(n);
  gsl_vector *k12 = gsl_vector_alloc(n);
  gsl_vector *temp_yt = gsl_vector_alloc(n);

  ftc(n,t,yt,k0,Tc,Tr,N);
  for (i = 0; i < n; i++) {
    gsl_vector_set(temp_yt,i,gsl_vector_get(yt,i) + gsl_vector_get(k0,i) * h / 2);
  }
  ftc(n,t+h/2,temp_yt,k12,Tc,Tr,N);
  for (i = 0; i < n; i++) {
    gsl_vector_set(yh,i,gsl_vector_get(yt,i) + gsl_vector_get(k12,i) * h);
  }
  for (i = 0; i < n; i++) {
    gsl_vector_set(err,i,(gsl_vector_get(k0,i) - gsl_vector_get(k12,i)) * h / 2);
  }
}

void driver (
	void f(int n, double t, gsl_vector *y, gsl_vector *dydx),
	double a,
	gsl_vector *ya,
	double b,
	double h,
	double acc,
	double eps
) {
  int n = ya->size;
  double x=a;
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector_memcpy(y,ya);
  gsl_vector *yh = gsl_vector_alloc(n);
  gsl_vector *err = gsl_vector_alloc(n);
  double s, normy, tol, error;

  while (x<b) {
    if (x+h>b) {
      h=b-x;
    }
    rkstep12(f,x,y,h,yh,err);
    s = 0;
    for (int i = 0; i < n; i++) {
      s+= gsl_vector_get(err,i)*gsl_vector_get(err,i);
    }
    error = sqrt(s);

    s=0;
    for (int i = 0; i < n; i++) {
      s+= gsl_vector_get(yh,i)*gsl_vector_get(yh,i);
    }
    normy=sqrt(s);

    FILE* stream = fopen("debug.txt","w");
    vector_print(stream,yh);
    fclose(stream);



    tol = (normy*eps+acc)*sqrt(h/(b-a));
    if (error<tol) {
      x+=h;
      y=yh;
      printf("%10g ",x);
      for (int i = 0; i < n; i++) {
        printf("%10g ",gsl_vector_get(yh,i));
      }
      printf("\n");
    }
    h*=pow(tol/error,0.25)*0.95;
  }
printf("\n\n\n\n");
}

void driverTC ( // TC driver
	void ftc(int n, double t, gsl_vector *y, gsl_vector *dydx,double Tc,double Tr, int N),
	double a,
	gsl_vector *ya,
	double b,
	double h,
	double acc,
	double eps,
  double Tc,
  double Tr,
  int N,
  FILE* stream2
) {
  int n = ya->size;
  double x=a;
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector_memcpy(y,ya);
  gsl_vector *yh = gsl_vector_alloc(n);
  gsl_vector *err = gsl_vector_alloc(n);
  double s, normy, tol, error;
  while (x<b) {
    if (x+h>b) {
      h=b-x;
    }
    rkstep12SIR(ftc,x,y,h,yh,err,Tc,Tr,N);
    s = 0;
    for (int i = 0; i < n; i++) {
      s+= gsl_vector_get(err,i)*gsl_vector_get(err,i);
    }
    error = sqrt(s);

    s=0;
    for (int i = 0; i < n; i++) {
      s+= gsl_vector_get(yh,i)*gsl_vector_get(yh,i);
    }
    normy=sqrt(s);

    tol = (normy*eps+acc)*sqrt(h/(b-a));
    if (error<tol) {
      x+=h;
      y=yh;
      fprintf(stream2,"%10g ",x);
      for (int i = 0; i < n; i++) {
        fprintf(stream2,"%10g ",gsl_vector_get(yh,i));
      }
      fprintf(stream2,"\n");
    }
    h*=pow(tol/error,0.25)*0.95;
  }
    fprintf(stream2,"\n\n\n\n");
}

int main() {

  // PART A.1 - A.3
  int m=2;
  gsl_vector *x0 = gsl_vector_alloc(m);
  gsl_vector_set(x0,0,0);
  gsl_vector_set(x0,1,1);

  void harm_osc (int n, double t, gsl_vector *x, gsl_vector *dx){
    gsl_vector_set(dx,1,-gsl_vector_get(x,0));
    gsl_vector_set(dx,0,gsl_vector_get(x,1));
  }

  double a=0.;
  double b=20.;

  double acc=0.0005;
  double eps=0.00005;
  double h=0.333;

  driver(&harm_osc,a,x0,b,h,acc,eps);

  // PART A.4

  void ftc(int n, double x, gsl_vector *y, gsl_vector *dydx, double Tc,double Tr, int N){
    gsl_vector_set(dydx,0,-1.0*(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc));
    assert(gsl_vector_get(dydx,0)<0);
    gsl_vector_set(dydx,1,(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-1.0*gsl_vector_get(y,1)/Tr);
    gsl_vector_set(dydx,2,(gsl_vector_get(y,1))/Tr);
  }

  int q = 3;
  a = 0.;
  b = 100.;
  double N=5.5e6;
  gsl_vector *yA = gsl_vector_alloc(q);
  gsl_vector_set(yA,1,661);
  gsl_vector_set(yA,2,5.5e5);
  gsl_vector_set(yA,0,N-gsl_vector_get(yA,2)-gsl_vector_get(yA,1));


  double Tc=5;
  double Tr=12;
  FILE* stream2 = fopen("out2.txt","w");
  driverTC(&ftc,a,yA,b,h,acc,eps,Tc,Tr,N,stream2);
  fclose(stream2);

  Tc=Tc/2;
  stream2 = fopen("out3.txt","w");
  driverTC(&ftc,a,yA,b,h,acc,eps,Tc,Tr,N,stream2);
  fclose(stream2);

  Tc=Tc/2;
  stream2 = fopen("out4.txt","w");
  driverTC(&ftc,a,yA,b,h,acc,eps,Tc,Tr,N,stream2);
  fclose(stream2);

  Tc=Tc/2;
  stream2 = fopen("out5.txt","w");
  driverTC(&ftc,a,yA,b,h,acc,eps,Tc,Tr,N,stream2);
  fclose(stream2);

  return 0;
}
