#include "stdio.h"
#include "math.h"

double my_exp(double x){
  if (x<0) {
    return 1/my_exp(-x);
  }
  if (x>1./8) {
    return pow(my_exp(x/2),2);
  }
  return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}


int main() {
  FILE* exp_out_stream = fopen("data.txt","w");
  double xmin = -5, xmax = 5;
  for (double x = xmin; x <=xmax; x += 1.0/8){
    fprintf(exp_out_stream,"%10g %10g %10g\n",x,exp(x),my_exp(x));
    }
  fclose(exp_out_stream);
  return 0;
}
