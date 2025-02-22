#ifndef NUMX_OPM_UNI_H
#define NUMX_OPM_UNI_H

#include <stdx/cap.h>

struct uni_bis_ops {
  double eps;
  double dta;

  struct vcap itr;
};

struct uni_gss_ops {
  struct vcap itr;
};

struct uni_fib_ops {
  double eps;

  struct vcap itr;
};

int uni_bis_min(double (*f)(double), double a, double b, double* r, struct uni_bis_ops o);
int uni_gss_min(double (*f)(double), double a, double b, double* r, struct uni_gss_ops o);
int uni_fib_min(double (*f)(double), double a, double b, double* r, struct uni_fib_ops o);

#endif  // NUMX_OPM_UNI_H