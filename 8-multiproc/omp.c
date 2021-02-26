#include<stdio.h>
#include<math.h>
#include<omp.h>

int main(){
	double x=0,y=-100,z=100;
	int N=(int)1e8;
#pragma omp parallel sections
	{
	#pragma omp section
		{
		for(int i=0;i<N;i++)x=cos(x);
		printf("x=%g cos(x)=%g\n",x,cos(x));
		}
	#pragma omp section
		{
		for(int i=0;i<N;i++)y=cos(y);
		printf("y=%g cos(y)=%g\n",y,cos(y));
		}
	#pragma omp section
		{
		for(int i=0;i<N;i++)z=cos(z);
		printf("z=%g cos(z)=%g\n",z,cos(z));
		}
	}
return 0;
}
