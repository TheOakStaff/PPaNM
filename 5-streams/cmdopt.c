#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>

int main (int argc, char *argv[]) {
	double epsilon=0.01; // default
	int npoints=9; // default
	while(1){
		int opt=getopt(argc,argv,"e:n:");
		if(opt==-1)break;
		switch(opt){
		case 'e': epsilon=atof(optarg); break;
		case 'n': npoints=atoi(optarg); break;
		default:
			fprintf(stderr,"usage: %s [ -e epsilon] {-n npoints]\n",argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	printf("epsilon=%g npoints=%i\n",epsilon,npoints);
exit(EXIT_SUCCESS);
return 0;
}
