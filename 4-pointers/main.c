#include<stdio.h>

void set0(double x){x=0; }

void set0p(double * x){ (*x)=0; }

int main(){
	double y=1;
	set0(y);
	printf("y after set0 = %g\n",y);
	set0p(&y);
	printf("y after set0p = %g\n",y);
	int n=5;
	double v[5];
	for(int i=0;i<n;i++){v[i]=i;} // i++: i=i+1
	int i=0; while (i<n) {
		printf("v[%d]=%g\n",i,v[i]);
		i++;
	}
return 0;
}
