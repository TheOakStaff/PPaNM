#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(){
	double x = 5;
	double r = tgamma(x);
	printf ("Gamma(%g) = %g\n",x,r);

	double y = 0.5;
	double t = j1(y);
	printf ("Bessel function of the first order of %g is equal to %g\n",y,t);

	double complex z = -2.0;
	double complex s = csqrt(z);
	printf ("The squareroot of -2 is %g + %gi\n",creal(s),cimag(s));

	double complex h = I * M_PI;
	double complex u = cpow(M_E, h);
	printf ("e to the power of %gi is equal to %f + %fi\n",cimag(h),creal(u),cimag(u));

	double complex j = cpow(M_E,I);
	printf ("e to the power of i is %f + %fi\n",creal(j),cimag(j));

	double complex k = cpow(I,M_E);
	printf ("i to the power of e is %f + %fi\n",creal(k),cimag(k));

	double complex m = cpow(I,I);
	printf ("i to the power of i is %f + %fi\n",creal(m),cimag(m));

	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;

	printf ("FL: %.25g\nDB: %.25lg\nLD: %.25Lg\n",x_float,x_double,x_long_double);

return 0;
}

