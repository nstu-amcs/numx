#include <assert.h>
#include <errno.h>
#include <numx/vec/iss.h>
#include <numx/vec/vec.h>

static int step(struct dmtx* m, struct vec* x, struct vec* f, struct vec* r, double omg) {
  int n = m->pps.n;
  int d = m->pps.d;

  double* fv = f->dat;
  double* xv = x->dat;
  double* rv = r->dat;

  int* la = m->la;
  double** ad = m->ad;

  for (int i = 0; i < n; ++i) {
    double pi = fv[i];

    for (int k = 0; k < d; ++k) {
      int j = la[k] + i;

      if (j < 0 || j >= n)
        continue;

      pi -= ad[i][k] * xv[j];
    }

    rv[i] = xv[i] + (omg / ad[i][0]) * pi;
  }

  return 0;
}

static inline void swap(struct vec* a, struct vec* b) {
  double* t = a->dat;

  a->dat = b->dat;
  b->dat = t;
}

int diss_jac_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_jac_ops o) {
  assert(m);
  assert(x);
  assert(f);

  assert(m->pps.n == x->n);
  assert(x->n == f->n);

  struct vec t;

  if (vec_new(&t, x->n))
    return -1;

  double nf = 0;
  double nt = 0;

  double res = 1;
  double eps = o.ops.eps;
  double omg = o.omg;

  int max = o.ops.max;

  vec_nrm(f, &nf);

  for (int k = 0; k < max && res >= eps; ++k) {
    step(m, x, f, &t, omg);
    swap(x, &t);

    mtx_vmlt(m, x, &t);
    vec_cmb(f, &t, &t, -1);
    vec_nrm(&t, &nt);

    res = nt / nf;

    if (o.ops.itr.call)
      o.ops.itr.call(o.ops.itr.ctx, 2, k, res);
  }

  vec_cls(&t);

  return 0;
}

int diss_rlx_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_rlx_ops o) {
  assert(m);
  assert(x);
  assert(f);

  assert(m->pps.n == x->n);
  assert(x->n == f->n);

  struct vec t;

  if (vec_new(&t, x->n))
    return -1;

  double nf = 0;
  double nt = 0;

  double res = 1;
  double eps = o.ops.eps;
  double omg = o.omg;

  int max = o.ops.max;

  vec_nrm(f, &nf);

  for (int k = 0; k < max && res >= eps; ++k) {
    step(m, x, f, x, omg);

    mtx_vmlt(m, x, &t);
    vec_cmb(f, &t, &t, -1);
    vec_nrm(&t, &nt);

    res = nt / nf;

    if (o.ops.itr.call)
      o.ops.itr.call(o.ops.itr.ctx, 2, k, res);
  }

  vec_cls(&t);

  return 0;
}

int diss_sor_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_sor_ops o) {
  assert(m);
  assert(x);
  assert(f);

  struct vec xo;
  double xon = 1;

  if (vec_new(&xo, x->n))
    return -1;

  double* fv = f->dat;
  double* xv = x->dat;

  double** ad = m->ad;
  int* la = m->la;

  for (int k = 1; k <= o.ops.max && xon > o.ops.eps; ++k) {
    vec_cpy(x, &xo);

    for (int i = 0; i < x->n; ++i) {
      double pi = fv[i];

      for (int e = 0; e < m->pps.d; ++e) {
        int j = la[e] + i;

        if (j < 0 || j >= x->n)
          continue;

        pi -= ad[i][e] * xv[j];
      }

      xv[i] = (1 - o.omg) * xv[i] + (o.omg / ad[i][0]) * pi;
    }

    vec_cmb(x, &xo, &xo, -1);
    vec_nrm(&xo, &xon);

    if (o.ops.itr.call)
      o.ops.itr.call(o.ops.itr.ctx, 2, k, xon);
  }

  return 0;
}