#include "pthread.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"


double len = 1e8;


int point(void* arg){
  int* seed = (int*)arg;
  printf("%d\n",seed );
  int in = 0;
  for (int i = 0; i < len; i++) {
    double x = (double)rand_r(seed) / (double)RAND_MAX;
    double y = (double)rand_r(seed) / (double)RAND_MAX;
    double l = sqrt(x*x + y*y);
    if (l <= 1) {
      in++;
    }
  }
  return(in);
}

int main() {
  int t = 1;
  int* pt = &t;
  int x = point(pt);
  printf("%d\n",x);
  double p = 4 * x / len;
  printf("Pi is roughly equal to %f\n",p);
  return 0;
}

/*
int point(void* arg){
  int* lst;
  int seed;
  int dummy;
  lst = (int*)arg;
  seed = lst[0];
  dummy = lst[1];
  printf("seed = %d\n",seed);
  printf("dummy = %d\n",dummy);
  int in = 0;
  for (int i = 0; i < len; i++) {
    double x = (double)rand_r(&seed) / (double)RAND_MAX;
    double y = (double)rand_r(&seed) / (double)RAND_MAX;
    double l = sqrt(x*x + y*y);
    if (l <= 1) {
      in++;
    }
  }
  return(in);
}

int main() {
  int t = 1;
  int* pt = &t;
  int a = 0;
  int* pa = &a;
  printf("%d\n",a);
  int* list[2];
  list[0] = pt;
  list[1] = pa;
  int x = point(list);
  printf("%d\n",x);
  double p = 4 * x / len;
  printf("Pi is roughly equal to %f\n",p);
  return 0;
}
*/
