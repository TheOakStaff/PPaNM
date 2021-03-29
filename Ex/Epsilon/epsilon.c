#include "stdio.h"
#include "math.h"
#include "limits.h"
#include "float.h"

int main() {

  printf("\nThe inbuild INT_MAX is = \t  %i\n", INT_MAX);

  int i = 1; while (i+1>i) {
    i++;
  }
  printf("While max int = \t\t  %i\n", i);

  i = 1;
  for (i = 1; i < i+1; i++) {
  }
  printf("For max int = \t\t\t  %i\n", i);

  i = 1;
  do {
    i++;
  } while(i+1>1);
  printf("Do max int = \t\t\t  %i\n", i);


  printf("\nThe inbuild INT_MIN is = \t %i\n", INT_MIN);

  i = 1; while (i-1<i) {
    i--;
  }
  printf("While min int = \t\t %i\n", i);

  i = 1;
  for (i = 1; i > i-1; i--) {
  }
  printf("For min int = \t\t\t %i\n", i);

  i = 1;
  do {
    i--;
  } while(i-1<1);
  printf("Do min int = \t\t\t %i\n", i);



  float wf = 1;
  while (1 + wf != 1) {
    wf/=2;
  }
  wf*=2;
  printf("\nFor Float_while_epsilon =\t  %g\n", wf);

  float ff = 1;
  for (ff = 1; 1 + ff != 1; ff/=2) {
  }
  ff*=2;
  printf("For Flaot_for_epsilon =\t\t  %g\n", ff);

  float df = 1;
  do {
    df/=2;
  } while(1 + df != 1);
  df*=2;
  printf("For Float_do_epsilon =\t\t  %g\n", df);


  double wd = 1;
  while (1 + wd != 1) {
    wd/=2;
  }
  wd*=2;
  printf("\nFor Double_while_epsilon =\t  %g\n", wd);

  double fd = 1;
  for (fd = 1; 1 + fd != 1; fd/=2) {
  }
  fd*=2;
  printf("For Double_for_epsilon =\t  %g\n", fd);

  double dd = 1;
  do {
    dd/=2;
  } while(1 + dd != 1);
  dd*=2;
  printf("For Double_do_epsilon =\t\t  %g\n", dd);

  long double wL = 1;
  while (1 + wL != 1) {
    wL/=2;
  }
  wL*=2;
  printf("\nFor Long_while_epsilon =\t  %Lg\n", wL);

  long double fL = 1;
  for (fL = 1; 1 + fL != 1; fL/=2) {
  }
  fL*=2;
  printf("For Long_for_epsilon =\t\t  %Lg\n", fL);

  long double dL = 1;
  do {
    dL/=2;
  } while(1 + dL != 1);
  dL*=2;
  printf("For Long_do_epsilon =\t\t  %Lg\n", dL);

  int max = INT_MAX/2;
  float suf = 0;
  float sdf = 0;
  for (int i = 1; i < max; i += 1) {
    suf += 1.0/i;
    sdf += 1.0/(max - i);
  }
  printf("The max is = %d\n", max);
  printf("The sum_up_flaot is = %f\n", suf);
  printf("The sum_down_flaot is = %f\n", sdf);

  printf("The difference is that when adding relatively equally sized number");
  printf("\ntogether the pressision is preserved\n");
  printf("When adding together numbers which are relatively far appart\n");
  printf("the pressision byond the larger number is lost\n");

  printf("\nIf we had infinite pressision the sums would converge as a\n");
  printf("functions of max\n\n");

  double sud = 0;
  double sdd = 0;
  for (int i = 1; i < max; i += 1) {
      sud += 1.0/i;
      sdd += 1.0/(max - i);
  }

  printf("The sum_up_double is = %g\n", sud);
  printf("The sum_down_double is = %g\n", sdd);
  printf("With higer pressision no data is lost");

  return 0;
}
