#include"komplex.h"
#include"stdio.h"

int main(){
  komplex z = {2,5};

  printf("First test: komplex_print\n");
  komplex_print("Z =",z);

/*
  printf("Second test: komplex_set\n");
  komplex* u_point = {0};
  komplex_set(u_point,3,6);
  komplex_print("U =",*u_point);
*/

  printf("Third test: komplex_new\n");
  komplex w = komplex_new(4,7);
  komplex_print("W =",w);


  printf("Fourth test: komplex_add\n");
  komplex zw = komplex_add(z,w);
  komplex_print("z + w =",zw);


  printf("Fifth test: komplex_sub\n");
  komplex wz = komplex_sub(w,z);
  komplex_print("w - z =",wz);

  return(0);
}
