#include <errno.h>
#include <numx/non/apx.h>
#include <numx/vec/dss.h>

// k a0 a1 a2 ... an
//   b0 b1 b2 ... bn
//   c0 c1 c2 ... cn
//   d0 d1 d2 ... dn
int apx_cub(struct vec* xv, struct imtx* km) {
  if (!xv || !km) {
    errno = EINVAL;
    return -1;
  }

  int n = xv->n - 1;
  int r = 0;

  struct imtx mm;
  struct vec gv;
  struct vec hv;

  if ((r = imtx_new(&mm, (struct imtx_pps){n + 1, n + 1})) && r)
    goto end;

  if ((r = vec_new(&gv, n + 1)) && r)
    goto end;

  if ((r = vec_new(&hv, n)) && r)
    goto end;

  double* g = gv.dat;
  double* h = hv.dat;
  double* x = xv->dat;

  double* a = km->dat[0];
  double* b = km->dat[1];
  double* c = km->dat[2];
  double* d = km->dat[3];

  for (int i = 0; i < n; ++i)
    h[i] = x[i + 1] - x[i];

  for (int i = 1; i < n; ++i)
    g[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

  mm.dat[0][0] = 1;
  mm.dat[n][n] = 1;

  for (int i = 1; i < n; ++i) {
    mm.dat[i][i - 1] = h[i - 1];
    mm.dat[i][i] = 2 * (h[i - 1] + h[i]);
    mm.dat[i][i + 1] = h[i];
  }

  struct vec cv = {
      .n = n + 1,
      .dat = c,
  };

  if ((r = dss_red_slv(&mm, &cv, &gv)) && r)
    goto end;

  for (int i = 0; i < n; ++i) {
    b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
    d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
  }

end:
  mtx_cls(&mm);
  vec_cls(&gv);
  vec_cls(&hv);

  return r;
}
