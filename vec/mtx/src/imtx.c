#include <errno.h>
#include <numx/vec/mtx.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int imtx_new(struct imtx* m, struct imtx_pps pps) {
  if (!m || pps.m == 0 || pps.n == 0) {
    errno = EINVAL;
    return -1;
  }

  m->pps = pps;
  m->dat = malloc(sizeof(double*) * pps.m);

  if (!m->dat) {
    errno = ENOMEM;
    return -1;
  }

  for (int i = 0; i < pps.m; ++i) {
    m->dat[i] = malloc(sizeof(double) * pps.n);

    if (!m->dat[i]) {
      for (int j = 0; j < i; ++j)
        free(m->dat[j]);

      free(m->dat);

      errno = ENOMEM;
      return -1;
    }

    memset(m->dat[i], 0, sizeof(double) * pps.n);
  }

  return 0;
}

int imtx_cls(struct imtx* m) {
  if (!m) {
    errno = EINVAL;
    return -1;
  }

  for (int i = 0; i < m->pps.m; ++i)
    free(m->dat[i]);

  free(m->dat);

  return 0;
}

int imtx_vmlt(struct imtx* m, struct vec* v, struct vec* r) {
  if (!m || !v || !r || m->pps.n != v->n || m->pps.m != r->n) {
    errno = EINVAL;
    return -1;
  }

  double** md = m->dat;
  double* vd = v->dat;
  double* rd = r->dat;

#pragma omp parallel for
  for (int i = 0; i < m->pps.m; ++i) {
    double s = 0;

#pragma omp parallel for reduction(+ : s)
    for (int j = 0; j < m->pps.n; ++j)
      s += md[i][j] * vd[j];

    rd[i] = s;
  }

  return 0;
}

int imtx_mmlt(struct imtx* a, struct imtx* b, struct imtx* r) {
  if (!a || !b || !r || a->pps.n != b->pps.m || a->pps.m != r->pps.m ||
      b->pps.n != r->pps.n) {
    errno = EINVAL;
    return -1;
  }

  double** ad = a->dat;
  double** bd = b->dat;
  double** rd = r->dat;

#pragma omp parallel for
  for (int i = 0; i < r->pps.m; ++i) {
#pragma omp parallel for
    for (int j = 0; j < r->pps.n; ++j) {
      double s = 0;

#pragma omp parallel for reduction(+ : s)
      for (int e = 0; e < a->pps.n; ++e)
        s += ad[i][e] * bd[e][i];

      rd[i][j] = s;
    }
  }

  return 0;
}
