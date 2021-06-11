#include "stdio.h"
#include "math.h"
#include "gsl/gsl_interp.h"

double binsearh (int n, double x[], double z){
  double start_index = 0;
  double end_index = n-1;
  while (start_index <= end_index){
     int middle = start_index + (end_index - start_index )/2;
     if (x[middle] < z && x[middle + 1] > z)
        return middle;
     if (x[middle] < z && x[middle + 1] <= z)
        start_index = middle + 1;
     else
        end_index = middle - 1;
  }
}

double linterp (int n,double x[],double y[], double z) {
  int i = binsearh(n,x,z);
  double slope = (y[i+1]-y[i])/(x[i+1]-x[i]);
  double value=y[i]+slope*(z-x[i]);
  return value;
}

double lininteg (int n, double x[],double y[], double z){
  double my_int = y[0];
  FILE* integral_out_stream = fopen("my_integral.txt","w");
  fprintf(integral_out_stream,"%d %f\n",0,0.);
    for (size_t i = 0; i < z; i++) {
      if (i + 1 < z) {
        double a = (y[i+1]-y[i])/(x[i+1]-x[i]);
        double b = y[i]-a*x[i];
        double F_upper = a/2.*(x[i+1]*x[i+1]) + b * x[i+1];
        double F_lower = a/2. * (x[i]*x[i]) + b * x[i];
        my_int = my_int + (F_upper - F_lower);
        fprintf(integral_out_stream,"%zu %f\n",i+1,my_int);
      }
      else {
        double a = (y[i+1]-y[i])/(x[i+1]-x[i]);
        double b = y[i]-a*x[i];
        double F_upper = a/2.*(z*z) + b * z;
        double F_lower = a/2. * (x[i]*x[i]) + b * x[i];
        my_int = my_int + (F_upper - F_lower);
        fprintf(integral_out_stream,"%f %f\n",z,my_int);
      }
    }
  fclose(integral_out_stream);
  return my_int;
}

int main() {
  int n = 5;
  double x[] = {0,1,2,3,4};
  double y[] = {0,4,-8,12,1};
  FILE* data1_out_stream = fopen("data.txt","w");
  for (size_t i = 0; i < n; i++) {
    fprintf(data1_out_stream,"%f,%f\n",x[i],y[i]);
  }
  fclose(data1_out_stream);
  double z = 4;
  double zy = linterp (n,x,y,z);
  double integ = lininteg(n,x,y,z);
  int i = binsearh(n,x,z);

// GSL PART

  gsl_interp* linear= gsl_interp_alloc(gsl_interp_linear,n);
  gsl_interp_init(linear,x,y,n);
  FILE* GSL_out_stream = fopen("GSL_data.txt","w");
  for (size_t i = 0; i < n; i++) {
    double interp_lin = gsl_interp_eval(linear,x,y,i,NULL);
    double integ_lin = gsl_interp_eval_integ(linear,x,y,x[0],i,NULL);
    fprintf(GSL_out_stream, "%zu %f %f\n",i,interp_lin,integ_lin);
  }
  fclose(GSL_out_stream);

  return 0;
}
