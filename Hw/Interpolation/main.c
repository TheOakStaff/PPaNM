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
  double a = (y[0+1]-y[0])/(x[0+1]-x[0]);
  double b = y[0]-a*x[0];
  double F_upper = a/2.*(z*z) + b * z;
  double F_lower = a/2. * (x[0]*x[0]) + b * x[0];
  double my_lin_integral = fabs(F_upper - F_lower);
  return my_lin_integral;
}

int main() {
  double x[] = {0,1,2,3,4};
  double y[] = {0,2,4,6,8};
  int n = 5;
  double z = 2.5;
  double zy = linterp (n,x,y,z);
  double integ = lininteg(n,x,y,z);
  int i = binsearh(n,x,z);
  gsl_interp *gsl_interpolation = gsl_interp_alloc(gsl_interp_linear,5);
  gsl_interp_init(gsl_interpolation,x,y,5);
  double gsl_eva = gsl_interp_eval(gsl_interpolation,x,y,z,gsl_interp_accel_alloc ());
  double gsl_int = gsl_interp_eval_integ(gsl_interpolation,x,y,0,z,gsl_interp_accel_alloc ());


  FILE* outtxt_out_stream = fopen("out.txt","w");
  fprintf(outtxt_out_stream,"%f is withing the interval of x[%d] to x[%d]\n\n",z,i,i+1);
  fprintf(outtxt_out_stream,"%f has the y value %f\n",z,zy);
  fprintf(outtxt_out_stream,"according to gsl %f has the y value %f\n\n",z,gsl_eva);
  fprintf(outtxt_out_stream,"The integral from %f to %f is equal to %f\n",x[0],z,integ);
  fprintf(outtxt_out_stream,"The gsl_integral from %f to %f is equal to %f\n",x[0],z,gsl_int);
  fclose(outtxt_out_stream);

  FILE* data_out_stream = fopen("my_data.txt","w");
  fprintf(data_out_stream,"%f,%f",z,zy);
  fclose(data_out_stream);
  return 0;
}
