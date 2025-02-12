#include <errno.h>
#include <numx/vec/iss.h>
#include <string.h>

static int siss_bcg_unc_slv(struct smtx* m, struct vec* x, struct vec* f, struct iss_bcg_opt o) {
  int n = m->pps.n;
  int c = 0;

  struct vec r;
  struct vec z;
  struct vec p;
  struct vec s;
  struct vec h;

  if (vec_new(&r, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&z, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&p, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&s, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&h, n)) {
    c = -1;
    goto end;
  }

  double alf = 0;
  double bet = 0;
  double omg = 0;

  double tmp = 0;
  double nrm = 0;

  mtx_vmlt(m, x, &r);
  vec_cmb(f, &r, &r, -1);

  vec_cpy(&r, &z);
  vec_cpy(&r, &p);

  for (int k = 1; k <= o.opt.max; ++k) {
    mtx_vmlt(m, &p, &h);

    vec_dot(&r, &z, &tmp);
    vec_dot(&h, &z, &bet);

    alf = tmp / bet;

    vec_cmb(&r, &h, &s, -alf);

    mtx_vmlt(m, &s, &r);
    vec_dot(&r, &s, &omg);
    vec_dot(&r, &r, &bet);

    omg = omg / bet;

    vec_cmb(x, &p, x, alf);
    vec_cmb(x, &s, x, omg);
    vec_cmb(&s, &r, &r, -omg);
    vec_nrm(&r, &nrm);

    if (o.opt.itr.call)
      o.opt.itr.call(o.opt.itr.ctx, 2, k, nrm);

    if (nrm < o.opt.eps)
      break;

    vec_dot(&r, &z, &bet);

    bet = (bet / tmp) * (alf / omg);

    vec_cmb(&p, &h, &p, -omg);
    vec_cmb(&r, &p, &p, bet);
  }

end:
  vec_cls(&r);
  vec_cls(&z);
  vec_cls(&p);
  vec_cls(&s);
  vec_cls(&h);

  return c;
}

static int siss_con_lslv(struct smtx* m, struct vec* x, struct vec* f) {
  int n = m->pps.n;

  int* mia = m->ia;
  int* mja = m->ja;

  double* mdr = m->dr;
  double* mlr = m->lr;

  double* xv = x->dat;
  double* fv = f->dat;

  for (int i = 0; i < n; ++i) {
    int lr0 = mia[i];
    int lr1 = mia[i + 1];

    double s = 0;

    for (int lr = lr0; lr < lr1; ++lr) {
      int j = mja[lr];

      s += mlr[lr] * xv[j];
    }

    xv[i] = (fv[i] - s) / mdr[i];
  }

  return 0;
}

static int siss_con_uslv(struct smtx* m, struct vec* x, struct vec* f) {
  int n = m->pps.n;

  int* mia = m->ia;
  int* mja = m->ja;

  double* mdr = m->dr;
  double* mur = m->ur;

  double* xv = x->dat;
  double* fv = f->dat;

  if (x != f)
    memcpy(xv, fv, sizeof(double) * n);

  for (int j = n - 1; j > -1; --j) {
    xv[j] = xv[j] / mdr[j];

    double x = xv[j];

    int ur0 = mia[j];
    int ur1 = mia[j + 1];

    for (int ur = ur0; ur < ur1; ++ur) {
      int ui = mja[ur];

      xv[ui] -= mur[ur] * x;
    }
  }

  return 0;
}

static int siss_bcg_con_slv(struct smtx* m, struct vec* x, struct vec* f, struct iss_bcg_opt o) {
  int n = m->pps.n;
  int c = 0;

  struct vec r;
  struct vec z;
  struct vec p;
  struct vec s;
  struct vec h;

  struct vec pt;
  struct vec st;

  if (vec_new(&r, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&z, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&p, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&s, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&h, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&pt, n)) {
    c = -1;
    goto end;
  }

  if (vec_new(&st, n)) {
    c = -1;
    goto end;
  }

  double alf = 0;
  double bet = 0;
  double omg = 0;

  double tmp = 0;
  double nrm = 0;

  mtx_vmlt(m, x, &r);
  vec_cmb(f, &r, &r, -1);

  vec_cpy(&r, &z);
  vec_cpy(&r, &p);

  for (int k = 1; k <= o.opt.max; ++k) {
    siss_con_lslv(o.con.sm, &pt, &p);
    siss_con_uslv(o.con.sm, &pt, &pt);

    mtx_vmlt(m, &pt, &h);

    vec_dot(&r, &z, &tmp);
    vec_dot(&h, &z, &bet);

    alf = tmp / bet;

    vec_cmb(&r, &h, &s, -alf);

    siss_con_lslv(o.con.sm, &st, &s);
    siss_con_uslv(o.con.sm, &st, &st);

    mtx_vmlt(m, &s, &r);
    vec_dot(&r, &s, &omg);
    vec_dot(&r, &r, &bet);

    omg = omg / bet;

    vec_cmb(x, &pt, x, alf);
    vec_cmb(x, &st, x, omg);
    vec_cmb(&s, &r, &r, -omg);
    vec_nrm(&r, &nrm);

    if (o.opt.itr.call)
      o.opt.itr.call(o.opt.itr.ctx, 2, k, nrm);

    if (nrm < o.opt.eps)
      break;

    vec_dot(&r, &z, &bet);

    bet = (bet / tmp) * (alf / omg);

    vec_cmb(&p, &h, &p, -omg);
    vec_cmb(&r, &p, &p, bet);
  }

end:
  vec_cls(&r);
  vec_cls(&z);
  vec_cls(&p);
  vec_cls(&s);
  vec_cls(&h);

  vec_cls(&st);
  vec_cls(&pt);

  return c;
}

int siss_bcg_slv(struct smtx* m, struct vec* x, struct vec* f, struct iss_bcg_opt o) {
  if (!m || !x || !f) {
    errno = EINVAL;
    return -1;
  }

  return o.con.sm ? siss_bcg_con_slv(m, x, f, o) : siss_bcg_unc_slv(m, x, f, o);
}
