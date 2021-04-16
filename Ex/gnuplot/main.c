#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double Erf(double);
double mygamma(double z);

int main(){
	FILE* erf_out_stream = fopen("data.txt","w");
	double xmin = -2, xmax = 2;
	for (double x = xmin; x <=xmax; x += 1.0/8){
		fprintf(erf_out_stream,"%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
		}
	fclose(erf_out_stream);

	FILE* gam_out_stream = fopen("gamdata.txt","w");
	double zmin = -5, zmax = 5;
//	fprintf(gam_out_stream, "test\n");
	for (double z = zmin; z <=zmax; z += 1.0/16){
			if (z <= 0 && floorf(z) == z) {
			}
			else {
				fprintf(gam_out_stream,"%10g %10g %10g %10g\n",z,tgamma(z),gsl_sf_gamma(z),mygamma(z));
			}
		}
	fclose(gam_out_stream);
return 0;
}
