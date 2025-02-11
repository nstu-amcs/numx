#include <errno.h>
#include <math.h>
#include <numx/vec/vec.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

int vec_new_ini(struct vec* v, int n, ...) {
  if (!v || n == 0) {
    errno = EINVAL;
    return -1;
  }

  v->n = n;
  v->dat = malloc(sizeof(double) * n);

  if (!v->dat) {
    errno = ENOMEM;
    return -1;
  }

  memset(v->dat, 0, sizeof(double) * n);

  va_list arg;
  va_start(arg, n);

  for (int i = 0; i < n; ++i) {
    double a = va_arg(arg, double);

    if (a == DBL_MAX)
      break;

    v->dat[i] = a;
  }

  va_end(arg);

  return 0;
}

int vec_cls(struct vec* v) {
  if (!v) {
    errno = EINVAL;
    return -1;
  }

  free(v->dat);

  return 0;
}

int vec_cmb(struct vec* a, struct vec* b, struct vec* r, double k) {
  if (!a || !b || !r || a->n != r->n || b->n != r->n) {
    errno = EINVAL;
    return -1;
  }

  int n = a->n;

  double* ad = a->dat;
  double* bd = b->dat;
  double* rd = r->dat;

#pragma omp parallel for
  for (int i = 0; i < n; ++i)
    rd[i] = ad[i] + k * bd[i];

  return 0;
}

int vec_dot(struct vec* a, struct vec* b, double* r) {
  if (!a || !b || !r) {
    errno = EINVAL;
    return -1;
  }

  int dim = a->n;

  double* ad = a->dat;
  double* bd = b->dat;

  double s = 0;

#pragma omp parallel for reduction(+ : s)
  for (int i = 0; i < dim; ++i)
    s += ad[i] * bd[i];

  *r = s;

  return 0;
}

int vec_nrm(struct vec* v, double* r) {
  if (!v || !r) {
    errno = EINVAL;
    return -1;
  }

  vec_dot(v, v, r);
  *r = sqrt(*r);

  return 0;
}

int vec_cpy(struct vec* s, struct vec* d) {
  if (!s || !d || s->n != d->n) {
    errno = EINVAL;
    return -1;
  }

  memcpy(d->dat, s->dat, sizeof(double) * s->n);

  return 0;
}

int vec_rst(struct vec* v) {
  if (!v) {
    errno = EINVAL;
    return -1;
  }

  memset(v->dat, 0, sizeof(double) * v->n);

  return 0;
}
