#include<stdio.h>
#include<math.h>
int main(){
	double x;
	int items;
	FILE* my_output_stream=fopen("outfile.txt","w");
	do{
		items = scanf("%lg",&x); // from stdin
		items=scanf("%lg",&x); // %lg = double precision, &x=into var(x)
		printf(stderr,"x=%g sin(x) = %g\n", x, sin(x)); // to stderr
		fprintf(my_output_stream,"x=%g sin(x) = %g\n", x, sin(x));
	}while(items!=EOF);
fclose(my_output_stream);
return 0;
}
