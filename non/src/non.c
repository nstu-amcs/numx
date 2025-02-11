#include <errno.h>
#include <math.h>
#include <numx/non/dif.h>
#include <numx/non/non.h>
#include <numx/vec/dss.h>
#include <numx/vec/mtx.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdx.h>

struct rec {
  int i;
  double v;
};

static int rec_cmp_dsc(void*, int n, ...) {
  if (n < 2)
    return 0;

  va_list arg;
  va_start(arg, n);

  struct rec* ar = (struct rec*)va_arg(arg, void*);
  struct rec* br = (struct rec*)va_arg(arg, void*);

  va_end(arg);

  if (fabs(ar->v) > fabs(br->v))
    return 1;

  if (fabs(ar->v) < fabs(br->v))
    return -1;

  return 0;
}

static int non_evo_exc(struct pcut* rec, struct pcut* fun, struct vec* xk, struct vec* fk, struct imtx* jk, struct non_new_opt* opt) {
  if (fun->len < xk->n) {
    errno = EINVAL;
    return -1;
  }

  struct rec** rp = (struct rec**)rec->dat;
  double (**fp)(struct vec*) = (double (**)(struct vec*))fun->dat;

  for (int i = 0; i < rec->len; ++i) {
    rp[i]->i = i;
    rp[i]->v = fp[i](xk);
  }

  if (rec->len > xk->n)
    if (cut_srt(rec))
      return -1;

  for (int i = 0; i < xk->n; ++i) {
    fk->dat[i] = -rp[i]->v;

    if (opt->jac)
      for (int j = 0; j < xk->n; ++j)
        jk->dat[i][j] = opt->jac->dat[rp[i]->i][j](xk);
    else
      for (int j = 0; j < xk->n; ++j)
        if (pdif(fp[rp[i]->i], j, opt->hop, xk, &jk->dat[i][j]))
          return -1;
  }

  return 0;
}

static int non_evo_con(struct pcut* rec, struct pcut* fun, struct vec* xk, struct vec* fk, struct imtx* jk, struct non_new_opt* opt) {
  if (non_evo_exc(rec, fun, xk, fk, jk, opt))
    return -1;

  if (rec->len == xk->n)
    return 0;

  struct rec** rp = (struct rec**)rec->dat;
  double (**fp)(struct vec*) = (double (**)(struct vec*))fun->dat;
  int l = xk->n - 1;

  if (rec->len > xk->n) {
    fk->dat[l] = 0;

    for (int j = 0; j < xk->n; ++j)
      jk->dat[l][j] = 0;

    for (int i = l; i < rec->len; ++i) {
      fk->dat[l] -= rp[i]->v * rp[i]->v;

      if (opt->jac)
        for (int j = 0; j < xk->n; ++j)
          jk->dat[l][j] += 2 * rp[i]->v * opt->jac->dat[rp[i]->i][j](xk);
      else {
        double pd = 0;

        for (int j = 0; j < xk->n; ++j) {
          if (pdif(fp[rp[i]->i], j, opt->hop, xk, &pd))
            return -1;

          jk->dat[l][j] += 2 * rp[i]->v * pd;
        }
      }
    }
  }

  return 0;
}

int non_new_slv(struct pcut* fun, struct vec* x, struct non_new_opt opt) {
  if (!fun || !x || fun->len < x->n) {
    errno = EINVAL;
    return -1;
  }

  struct pcut rec;

  if (cut_new(&rec))
    goto err;

  if (cut_exp(&rec, fun->len))
    goto err;

  rec.ctl = true;
  rec.cmp.call = &rec_cmp_dsc;

  for (int i = 0; i < fun->len; ++i) {
    struct rec* r = malloc(sizeof(struct rec));

    if (!r) {
      errno = ENOMEM;
      goto err;
    }

    if (cut_add(&rec, r))
      goto err;
  }

  struct vec fk;
  struct vec dk;
  struct imtx jk;

  double nrm = 0;
  int dim = x->n;

  if (vec_new(&fk, dim))
    goto err;

  if (vec_new(&dk, dim))
    goto err;

  struct imtx_pps p = {dim, dim};

  if (mtx_new(&jk, p))
    goto err;

  switch (opt.mod) {
    case EXC:
      if (non_evo_exc(&rec, fun, x, &fk, &jk, &opt))
        goto err;

      break;
    case CON:
      if (non_evo_con(&rec, fun, x, &fk, &jk, &opt))
        goto err;

      break;
  }

  if (opt.itr) {
    opt.itr->k = 0;
    opt.itr->x = x;
    opt.itr->del = -1;

    if (vec_nrm(&fk, &opt.itr->err))
      goto err;
  }

  if (opt.cbk)
    opt.cbk->call(opt.cbk->ctx, 1, opt.itr);

  for (int k = 1; k <= opt.hem; ++k) {
    if (dss_red_slv(&jk, &dk, &fk))
      goto err;

    if (vec_cmb(x, &dk, x, 1))
      goto err;

    if (vec_nrm(&dk, &nrm))
      goto err;

    switch (opt.mod) {
      case EXC:
        if (non_evo_exc(&rec, fun, x, &fk, &jk, &opt))
          goto err;

        break;
      case CON:
        if (non_evo_con(&rec, fun, x, &fk, &jk, &opt))
          goto err;

        break;
    }

    if (opt.itr) {
      opt.itr->k = k;
      opt.itr->del = nrm;

      if (vec_nrm(&fk, &opt.itr->err))
        goto err;
    }

    if (opt.cbk)
      opt.cbk->call(opt.cbk->ctx, 1, opt.itr);

    if (nrm < opt.eps)
      break;
  }

  int r = 0;

  goto end;

err:
  r = -1;

end:
  cut_cls(&rec);
  vec_cls(&fk);
  vec_cls(&dk);
  mtx_cls(&jk);

  return r;
}
