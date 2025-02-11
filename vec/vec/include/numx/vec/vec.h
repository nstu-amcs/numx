#ifndef NUMX_VEC_VEC_H
#define NUMX_VEC_VEC_H

#include <float.h>

typedef unsigned int uint;

struct vec {
  int n;

  double* dat;
};

int vec_new_ini(struct vec* v, int n, ...);
int vec_cls(struct vec* v);

int vec_cmb(struct vec* a, struct vec* b, struct vec* r, double k);
int vec_crs(struct vec* a, struct vec* b, struct vec* r);
int vec_dot(struct vec* a, struct vec* b, double* r);
int vec_nrm(struct vec* v, double* r);

int vec_cpy(struct vec* s, struct vec* d);
int vec_rst(struct vec* v);

#define vec_new(v, ...) vec_new_ini(v, __VA_ARGS__, DBL_MAX)

#endif  // NUMX_VEC_VEC_H
