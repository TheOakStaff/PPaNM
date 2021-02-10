#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
	int n=10;
	double x,y=1;
	complex z;
	z = csqrt(-2);
	printf("%f%+fi\n",creal(z),cimag(z));
}
