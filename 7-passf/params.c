#include<stdio.h>
#include<math.h>

double k=1;

double g(double x){ return sin(k*x); }

typedef struct {double (*f)(double,void*); void* params;} myfun;

double f(double x,void* params){
	double k=*(double*)params;
	return sin(k*x);
}

void print_table(double f(double),double a, double b, double dx){
	for(double x=a;x<=b;x+=dx)
		printf("%10g %10g\n",x,f(x));
}

void print_table_params(myfun F,double a, double b, double dx){
	for(double x=a;x<=b;x+=dx)
		printf("%10g %10g\n",x,F.f(x,F.params));
}

int main(){
	double w=1;
	myfun F;
	F.f=&f;
	F.params=(void*)&w;
	print_table_params(F,0,M_PI,M_PI/8);
return 0;
}
