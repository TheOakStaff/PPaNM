#include "pthread.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

double len = 1e8;

int carlo(void* arg){
  int s = rand();
  // calling rand inside of a thread is safe since 
  // only the thread in question has acces to the rand value.
  int* seed = &s;
  int* in = (int*)arg;
  for (int i = 0; i < len; i++) {
    double x = (double)rand_r(seed) / (double)RAND_MAX;
    double y = (double)rand_r(seed) / (double)RAND_MAX;
    double l = sqrt(x*x + y*y);
    if (l <= 1) {
      *in = *in + 1;
    }
  }
  return(&in);
}

int main() {
  int t1 = 0, t2=0,t3=0,t4=0;
  pthread_t thread1, thread2, thread3, thread4;
	pthread_attr_t* attributes = NULL;
	pthread_create( &thread1, attributes, carlo, (void*)&t1 );
	pthread_create( &thread2, attributes, carlo, (void*)&t2 );
	pthread_create( &thread3, attributes, carlo, (void*)&t3 );
	pthread_create( &thread4, attributes, carlo, (void*)&t4 );
  void* returnvalue = NULL;
	pthread_join(thread1,returnvalue);
	pthread_join(thread2,returnvalue);
	pthread_join(thread3,returnvalue);
	pthread_join(thread4,returnvalue);
  printf("1 gives = %d\n2 gives = %d\n3 gives = %d\n4 gives = %d\n",t1,t2,t3,t4 );
  int tot = t1 + t2 + t3 + t4;
  double p = 4 * (double)tot / (double)(4*len);
  printf("Pi is roughly equal to %f\n",p);
  return 0;
}
