#include "stdio.h"
#include "komplex.h"

void komplex_print(char* s, komplex z) {
  printf("%s (%g,%g)\n", s, z.re, z.im);
}

void komplex_set(komplex* z, double x, double y) {
  (*z).re = x;
  (*z).im = y;
}

komplex komplex_new(double x, double y) {
  komplex z = {x, y};
  return z;
}

komplex komplex_add(komplex x, komplex y) {
  komplex add_result = {x.re + y.re, x.im + y.im};
  return add_result;
}

komplex komplex_sub(komplex x, komplex y) {
  komplex sub_result = {x.re - y.re, x.im - y.im};
  return sub_result;
}
