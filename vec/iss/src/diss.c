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

int diss_jac_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_jac_opt o) {
  if (!m || !x || !f) {
    errno = EINVAL;
    return -1;
  }

  struct vec t;

  if (vec_new(&t, x->n))
    return -1;

  double nf = 0;
  double nt = 0;

  double res = 1;
  double eps = o.opt.eps;
  double omg = o.omg;

  int max = o.opt.max;

  vec_nrm(f, &nf);

  for (int k = 0; k < max && res >= eps; ++k) {
    step(m, x, f, &t, omg);
    swap(x, &t);

    mtx_vmlt(m, x, &t);
    vec_cmb(f, &t, &t, -1);
    vec_nrm(&t, &nt);

    res = nt / nf;

    if (o.opt.itr.call)
      o.opt.itr.call(o.opt.itr.ctx, 2, k, res);
  }

  vec_cls(&t);

  return 0;
}

int diss_rlx_slv(struct dmtx* m, struct vec* x, struct vec* f, struct iss_rlx_opt o) {
  if (!m || !x || !f) {
    errno = EINVAL;
    return -1;
  }

  struct vec t;

  if (vec_new(&t, x->n))
    return -1;

  double nf = 0;
  double nt = 0;

  double res = 1;
  double eps = o.opt.eps;
  double omg = o.omg;

  int max = o.opt.max;

  vec_nrm(f, &nf);

  for (int k = 0; k < max && res >= eps; ++k) {
    step(m, x, f, x, omg);

    mtx_vmlt(m, x, &t);
    vec_cmb(f, &t, &t, -1);
    vec_nrm(&t, &nt);

    res = nt / nf;

    if (o.opt.itr.call)
      o.opt.itr.call(o.opt.itr.ctx, 2, k, res);
  }

  vec_cls(&t);

  return 0;
}