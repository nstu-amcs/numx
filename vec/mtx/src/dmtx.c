#include <errno.h>
#include <numx/vec/mtx.h>
#include <stdlib.h>
#include <string.h>

int dmtx_new(struct dmtx* m, struct dmtx_pps pps) {
  if (!m || pps.n == 0 || pps.d == 0) {
    errno = EINVAL;
    return -1;
  }

  m->pps = pps;
  m->la = malloc(sizeof(int) * pps.d);

  if (!m->la) {
    errno = ENOMEM;
    return -1;
  }

  memset(m->la, 0, sizeof(int) * pps.d);
  m->ad = malloc(sizeof(double*) * pps.n);

  if (!m->ad) {
    free(m->la);

    errno = ENOMEM;
    return -1;
  }

  for (int i = 0; i < pps.n; ++i) {
    m->ad[i] = malloc(sizeof(double) * pps.d);

    if (!m->ad[i]) {
      for (int j = 0; j < i; ++j)
        free(m->ad[j]);

      free(m->ad);
      free(m->la);

      errno = ENOMEM;
      return -1;
    }

    memset(m->ad[i], 0, sizeof(double) * pps.d);
  }

  return 0;
}

int dmtx_cls(struct dmtx* m) {
  if (!m) {
    errno = EINVAL;
    return -1;
  }

  for (int i = 0; i < m->pps.n; ++i)
    free(m->ad[i]);

  free(m->ad);
  free(m->la);

  return 0;
}

int dmtx_vmlt(struct dmtx* m, struct vec* x, struct vec* f) {
  if (!m || !x || !f) {
    errno = EINVAL;
    return -1;
  }

  int n = m->pps.n;
  int d = m->pps.d;

  double* xv = x->dat;
  double* fv = f->dat;

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    double* rv = m->ad[i];

    for (int e = 0; e < d; ++e) {
      int j = i + m->la[e];

      if (j < 0)
        continue;

      fv[i] += rv[e] * xv[j];
    }
  }

  return 0;
}