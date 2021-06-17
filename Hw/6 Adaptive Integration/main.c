#include "stdio.h"
#include "math.h"
#include "assert.h"


double adapt24 (double f (double ), double a, double b,
 double acc, double eps, double f2, double f3, int nrec){
	assert(nrec<1000000);
	double f1 = f(a+(b-a)/6);
  double f4 = f(a+5*(b-a)/6);
	double Q=(2*f1 + f2 + f3 + 2*f4)/6*(b - a);
  double q = (f1 + f4 + f2 + f3)/4*(b - a);
	double tolerance = acc + eps*fabs(Q);
  double error = fabs(Q-q) ;
	if (error < tolerance) return Q;
	else {
		double Q1=adapt24(f, a, (a+b)/2,acc/sqrt(2.), eps, f1, f2, nrec+1);
		double Q2=adapt24(f, (a+b)/2, b, acc/sqrt(2.), eps, f3, f4, nrec+1);
		return Q1+Q2;
  }
}
double adapt (double f(double), double a, double b, double acc, double eps){
	double f2=f(a+2*(b-a)/6);
  double f3=f(a+4*(b-a)/6);
	int nrec=0;
	return adapt24(f ,a ,b ,acc ,eps ,f2 ,f3, nrec);
}
int main() {
	int calls=0; double a=0, b=1, acc=0.0001, eps = 0.00001;
	double f(double x){
    calls++;
    return sqrt(x)-2./3.;
  };
	double Q = adapt(f, a ,b ,acc ,eps);
	printf("\nIntegral of sqrt(x)-2/3 = %g\nTo achive sufficent accuracy the fucntion was called = %d\n", Q, calls);
  printf("The true integral of sqrt(x)-2/3 = 0\n\n");

  calls=0;
  double g(double x){
    calls++;
    return 4*sqrt(1-pow(x,2)) - M_PI;
  };
  double P = adapt(g, a ,b ,acc ,eps);
  printf("Integral of 4*sqrt(1-x²) - π = %g\nTo achive sufficent accuracy the fucntion was called = %d\n", P, calls);
  printf("The true integral of 4*sqrt(1-x²) - π = 0\n\n");

  calls=0;
  double h(double x){
    calls++;
    return pow(x,2)- pow(x,3);
  };
  double H = adapt(h, a ,b ,acc ,eps);
  printf("Integral of x²-x³ = %g\nTo achive sufficent accuracy the fucntion was called = %d\n", H, calls);
  printf("The true integral of x²-x³ = 0.83(3)\n");
return 0; }
