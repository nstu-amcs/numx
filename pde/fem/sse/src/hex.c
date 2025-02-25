#include "const.h"

#include <fem/sse/hex.h>
#include <vec/mtx.h>

#include <errno.h>
#include <stdlib.h>

int hex_ini(struct hex** h) {
  if (!h) {
    errno = EINVAL;
    return -1;
  }

  struct hex* hex = malloc(sizeof(struct hex));

  if (!hex) {
    errno = ENOMEM;
    return -1;
  }

  hex->vtx[0] = 0;
  hex->vtx[1] = 0;
  hex->vtx[2] = 0;
  hex->vtx[3] = 0;
  hex->vtx[4] = 0;
  hex->vtx[5] = 0;
  hex->vtx[6] = 0;
  hex->vtx[7] = 0;

  hex->loc.b = 0;
  hex->loc.m = 0;

  hex->pps.gam = 0;
  hex->pps.lam = 0;

  hex->pps.f = 0;

  *h = hex;

  return 0;
}

int hex_cls(struct hex** h) {
  if (!h || !(*h)) {
    errno = EINVAL;
    return -1;
  }

  struct hex* hex = *h;

  if (hex->loc.m)
    mtx_cls(hex->loc.m);

  if (hex->loc.b)
    vec_cls(hex->loc.b);

  *h = 0;
  free(hex);

  return 0;
}

int hex_new(struct hex* h) {
  if (!h) {
    errno = EINVAL;
    return -1;
  }

  h->loc.m = mtx_new(8);
  h->loc.b = vec_new(8);

  return 0;
}

int hex_sget(struct hex* h, const char* buf, double (**fun)(struct vtx*)) {
  if (!h || !buf || !fun || buf[0] != 'h') {
    errno = EINVAL;
    return -1;
  }

  buf += 1;

  int pos, f;

  for (int j = 0, v; j < 8; ++j) {
    sscanf(buf, "%d%n", &v, &pos);
    buf += pos;

    h->vtx[j] = v - 1;
  }

  sscanf(buf, " | %d | %lf | %lf%n", &f, &h->pps.lam, &h->pps.gam, &pos);
  buf += pos;

  h->pps.f = fun[f - 1];

  return 0;
}

int hex_evo(struct hex* h, struct vtx** v) {
  double lam = h->pps.lam;
  double gam = h->pps.gam;

  double hx = v[h->vtx[1]]->x - v[h->vtx[0]]->x;
  double hy = v[h->vtx[2]]->y - v[h->vtx[0]]->y;
  double hz = v[h->vtx[4]]->z - v[h->vtx[0]]->z;

  double gx[2][2];
  double gy[2][2];
  double gz[2][2];

  double mx[2][2];
  double my[2][2];
  double mz[2][2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      gx[i][j] = G[i][j] / hx;
      gy[i][j] = G[i][j] / hy;
      gz[i][j] = G[i][j] / hz;

      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;
    }

  double** mloc = h->loc.m->v;
  double* bloc = h->loc.b->vp;
  double dec[8];

  //  static int sign[2] = {-1, 1};
  //
  //  for (int i = 0; i < 8; ++i)
  //    dec[i] = h->pps.f(v[h->vtx[i]]) * sign[MU[i]] / hx;

  for (int i = 0; i < 8; ++i)
    dec[i] = h->pps.f(v[h->vtx[i]]);

  for (int i = 0; i < 8; ++i) {
    bloc[i] = 0;

    int mui = MU[i];
    int nui = NU[i];
    int tti = TT[i];

    for (int j = 0; j < 8; ++j) {
      int muj = MU[j];
      int nuj = NU[j];
      int ttj = TT[j];

      double mxl = mx[mui][muj];
      double myl = my[nui][nuj];
      double mzl = mz[tti][ttj];

      double gxl = gx[mui][muj];
      double gyl = gy[nui][nuj];
      double gzl = gz[tti][ttj];

      mloc[i][j] = gam * mxl * myl * mzl;
      mloc[i][j] += lam * (gxl * myl * mzl + mxl * gyl * mzl + mxl * myl * gzl);

      bloc[i] += dec[j] * mxl * myl * mzl;
      // bloc[i] += dec[j] * (hx / 2) * myl * mzl;
    }
  }

  return 0;
}

int hex_mov(struct hex* h, struct mtx_csj* a, struct vec* b) {
  double* bloc = h->loc.b->vp;
  double** mloc = h->loc.m->v;

  for (int i = 0; i < 8; ++i) {
    int iglob = h->vtx[i];

    b->vp[iglob] += bloc[i];

    for (int j = 0; j < 8; ++j) {
      int jglob = h->vtx[j];

      if (jglob == iglob) {
        a->dr[iglob] += mloc[i][j];
        continue;
      }

      if (iglob > jglob) {
        int lr0 = a->ia[iglob];
        int lr1 = a->ia[iglob + 1];

        for (int lr = lr0; lr < lr1; ++lr) {
          int lj = a->ja[lr];

          if (lj == jglob) {
            a->lr[lr] += mloc[i][j];
            break;
          }
        }
      } else {
        int ur0 = a->ia[jglob];
        int ur1 = a->ia[jglob + 1];

        for (int ur = ur0; ur < ur1; ++ur) {
          int ui = a->ja[ur];

          if (ui == iglob) {
            a->ur[ur] += mloc[i][j];
            break;
          }
        }
      }
    }
  }

  return 0;
}
