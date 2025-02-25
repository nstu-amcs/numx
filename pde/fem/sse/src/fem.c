#include "const.h"

#include <fem/sse/fem.h>
#include <vec/iss_csj.h>
#include <ext/sll.h>

#include <errno.h>
#include <stdlib.h>
#include <string.h>

const double G[2][2] = {{1.0, -1.0}, {-1.0, 1.0}};
const double M[2][2] = {{2.0 / 6.0, 1.0 / 6.0}, {1.0 / 6.0, 2.0 / 6.0}};
const double X[2][2] = {{-0.5, 0.5}, {-0.5, 0.5}};

const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

const double C = 10e10;

int fem_ini(struct fem** h) {
  if (!h) {
    errno = EINVAL;
    return -1;
  }

  struct fem* fem = malloc(sizeof(struct fem));

  if (!fem) {
    errno = ENOMEM;
    return -1;
  }

  fem->vs = 0;
  fem->hs = 0;
  fem->fs = 0;

  fem->a = 0;
  fem->b = 0;

  fem->vtx = 0;
  fem->hex = 0;
  fem->fce = 0;

  fem->pps.fun = 0;

  *h = fem;

  return 0;
}

int fem_cls(struct fem** h) {
  if (!h || !(*h)) {
    errno = EINVAL;
    return -1;
  }

  struct fem* fem = *h;

  if (fem->a)
    mtx_csj_cls(fem->a);

  if (fem->b)
    vec_cls(fem->b);

  if (fem->vtx) {
    for (int i = 0; i < fem->vs; ++i)
      vtx_cls(&fem->vtx[i]);

    free(fem->vtx);
  }

  if (fem->hex) {
    for (int i = 0; i < fem->hs; ++i)
      hex_cls(&fem->hex[i]);

    free(fem->hex);
  }

  if (fem->fce) {
    for (int i = 0; i < fem->fs; ++i)
      fce_cls(&fem->fce[i]);

    free(fem->fce);
  }

  free(fem);
  *h = 0;

  return 0;
}

int fem_fget(struct fem* fem, FILE* f) {
  char buf[0xff];

  fgets(buf, sizeof(buf), f);
  sscanf(buf, "#vertex [%d]", &fem->vs);

  fem->vtx = malloc(sizeof(struct vtx*) * fem->vs);

  for (int i = 0; i < fem->vs; ++i) {
    fgets(buf, sizeof(buf), f);

    vtx_ini(&fem->vtx[i]);
    vtx_sget(fem->vtx[i], buf);
  }

  fgets(buf, sizeof(buf), f);
  fgets(buf, sizeof(buf), f);
  sscanf(buf, "#hexahedron [%d]", &fem->hs);

  fem->hex = malloc(sizeof(struct hex*) * fem->hs);

  for (int i = 0; i < fem->hs; ++i) {
    fgets(buf, sizeof(buf), f);

    hex_ini(&fem->hex[i]);
    hex_new(fem->hex[i]);
    hex_sget(fem->hex[i], buf, fem->pps.fun);
  }

  fgets(buf, sizeof(buf), f);
  fgets(buf, sizeof(buf), f);
  sscanf(buf, "#face [%d]", &fem->fs);

  fem->fce = malloc(sizeof(struct fce*) * fem->fs);

  for (int i = 0; i < fem->fs; ++i) {
    fgets(buf, sizeof(buf), f);

    fce_ini(&fem->fce[i]);
    fce_sget(fem->fce[i], buf, fem->pps.fun);
  }

  return 0;
}

static int cmp(int a, int b) {
  if (a < b)
    return 1;

  if (a > b)
    return -1;

  return 0;
}

int fem_evo(struct fem* fem) {
  int ne = 0;

  struct isll** l = malloc(sizeof(struct isll*) * fem->vs);

  for (int i = 0; i < fem->vs; ++i) {
    sll_ini(&l[i]);
    sll_cmp(l[i], &cmp);
    sll_dup(l[i], false);
  }

  for (int i = 0; i < fem->hs; ++i) {
    struct hex* h = fem->hex[i];

    hex_evo(h, fem->vtx);

    for (int j = 0; j < 8; ++j) {
      int a = h->vtx[j];

      for (int k = 0; k < 8; ++k) {
        int b = h->vtx[k];

        if (b < a && !sll_add(l[a], b))
          ne += 1;
      }
    }
  }

  for (int i = 0; i < fem->fs; ++i)
    fce_evo(fem->fce[i], fem->vtx);

  fem->a = mtx_csj_new(fem->vs, ne);
  fem->b = vec_new(fem->vs);

  for (int i = 0, e = 0; i < fem->vs; ++i) {
    fem->a->ia[i] = e;

    struct isln* j = 0;

    while (sll_next(l[i], &j)) {
      fem->a->ja[e] = j->e;
      e = e + 1;
    }
  }

  fem->a->ia[fem->vs] = ne;

  for (int i = 0; i < fem->vs; ++i)
    sll_cls(&l[i]);

  free(l);

  return 0;
}

int fem_asm(struct fem* fem) {
  struct isll* dir;

  sll_ini(&dir);

  for (int i = 0; i < fem->hs; ++i)
    hex_mov(fem->hex[i], fem->a, fem->b);

  for (int i = 0; i < fem->fs; ++i) {
    struct fce* f = fem->fce[i];

    switch (f->cnd.type) {
      case DIR:
        sll_add(dir, i);
        continue;
      case NEU:
      case ROB:
        fce_mov(f, fem->a, fem->b);
        continue;
    }
  }

  struct isln* i = 0;

  while (sll_next(dir, &i)) {
    struct fce* f = fem->fce[i->e];

    for (int j = 0; j < 4; ++j) {
      struct vtx* v = fem->vtx[f->vtx[j]];

      fem->a->dr[f->vtx[j]] = C;
      fem->b->vp[f->vtx[j]] = C * f->cnd.pps.dir.tmp(v);
    }
  }

  sll_cls(&dir);

  return 0;
}

int fem_slv(struct fem* fem, struct vec* q) {
  struct iss_pps pps = {1e-5, 10000};
  struct iss_res res = {0, 0};

  iss_csj_bcg_slv(fem->a, q, fem->b, &pps, &res, NULL);

  return 0;
}

int fem_get(struct fem* fem, struct vec* q, struct vtx* v, uint n, double* r) {
  struct hex* h = fem->hex[n];

  double xf = fem->vtx[h->vtx[0]]->x;
  double xs = fem->vtx[h->vtx[1]]->x;
  double xh = xs - xf;

  double yf = fem->vtx[h->vtx[0]]->y;
  double ys = fem->vtx[h->vtx[2]]->y;
  double yh = ys - yf;

  double zf = fem->vtx[h->vtx[0]]->z;
  double zs = fem->vtx[h->vtx[4]]->z;
  double zh = zs - zf;

  double s = 0;

  for (int i = 0; i < 8; ++i) {
    double p = 1;

    p *= MU[i] ? (v->x - xf) / xh : (xs - v->x) / xh;
    p *= NU[i] ? (v->y - yf) / yh : (ys - v->y) / yh;
    p *= TT[i] ? (v->z - zf) / zh : (zs - v->z) / zh;

    s += p * q->vp[h->vtx[i]];
  }

  *r = s;

  return 0;
}
