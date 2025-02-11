#include <errno.h>
#include <numx/vec/mtx.h>
#include <stdlib.h>
#include <string.h>

int jmtx_new(struct jmtx* m, struct jmtx_pps pps) {
  if (!m || pps.n == 0 || pps.m == 0) {
    errno = EINVAL;
    return -1;
  }

  m->pps = pps;
  m->dat = malloc(sizeof(double (*)(struct vec*)) * pps.m);

  if (!m->dat) {
    errno = ENOMEM;
    return -1;
  }

  for (int i = 0; i < pps.m; ++i) {
    m->dat[i] = malloc(sizeof(double (*)(struct vec*)) * pps.n);

    if (!m->dat[i]) {
      for (int j = 0; j < i; ++j)
        free(m->dat[j]);

      free(m->dat);

      errno = ENOMEM;
      return -1;
    }

    memset(m->dat[i], 0, sizeof(double (*)(struct vec*)) * pps.n);
  }

  return 0;
}

int jmtx_cls(struct jmtx* m) {
  if (!m) {
    errno = EINVAL;
    return -1;
  }

  for (int i = 0; i < m->pps.m; ++i)
    free(m->dat[i]);

  free(m->dat);

  return 0;
}
