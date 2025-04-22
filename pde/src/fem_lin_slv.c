#include "fem_lin.h"

#define HOP 0.0001

double fem_lin_apx(struct apx_fun_ctx *ctx, struct vec *vtx) {
  assert(ctx);
  assert(vtx);

  struct vec *v = ctx->sim->msh->vtx.dat;
  struct hxd *h = &ctx->sim->msh->hxd.dat[ctx->hxd];

  double *w = ctx->wgt->dat;

  int v0 = h->vtx[0];
  int v7 = h->vtx[7];

  double x1 = v[v0].dat[0];
  double x2 = v[v7].dat[0];
  double y1 = v[v0].dat[1];
  double y2 = v[v7].dat[1];
  double z1 = v[v0].dat[2];
  double z2 = v[v7].dat[2];

  double hm = (x2 - x1) * (y2 - y1) * (z2 - z1);

  double x = vtx->dat[0];
  double y = vtx->dat[1];
  double z = vtx->dat[2];
  double r = 0;

  r += w[h->vtx[0]] * (x2 - x) * (y2 - y) * (z2 - z) / hm;
  r += w[h->vtx[1]] * (x - x1) * (y2 - y) * (z2 - z) / hm;
  r += w[h->vtx[2]] * (x2 - x) * (y - y1) * (z2 - z) / hm;
  r += w[h->vtx[3]] * (x - x1) * (y - y1) * (z2 - z) / hm;
  r += w[h->vtx[4]] * (x2 - x) * (y2 - y) * (z - z1) / hm;
  r += w[h->vtx[5]] * (x - x1) * (y2 - y) * (z - z1) / hm;
  r += w[h->vtx[6]] * (x2 - x) * (y - y1) * (z - z1) / hm;
  r += w[h->vtx[7]] * (x - x1) * (y - y1) * (z - z1) / hm;

  return r;
}

static int ell_slv(struct sim *sim);
static int pbc_slv(struct sim *sim);

int fem_lin_slv(struct sim *sim) {
  assert(sim);

  sim->slv->apx = fem_lin_apx;

  switch (sim->mod) {
  case SIM_ELL:
    return ell_slv(sim);
  case SIM_PBC:
    return pbc_slv(sim);
  case SIM_HYP:
    return pbc_slv(sim);
  }

  return 0;
}

static int sys_slv(struct sim *sim);

static int ell_slv(struct sim *sim) {
  int r = 0;

  struct fem *fem = (struct fem *)sim->slv;
  struct vec wgt;

  if ((r = vec_new(&wgt, fem->prv.vec.n)))
    goto end;

  fem->slv.run.ti = 0;
  fem->slv.run.tv = 0;
  fem->slv.run.tb = 1;
  fem->slv.run.wgt[0] = &wgt;

  if ((r = sys_slv(sim)))
    goto end;

  if (fem->slv.itr.run)
    fem->slv.itr.run(fem->slv.itr.ctx, sim);

  if (sim->ops.exp.put)
    sim->ops.exp.put(sim);

end:
  vec_cls(&wgt);

  return r;
}

static int pbc_i1s_slv(struct sim *sim);

static int pbc_slv(struct sim *sim) {
  int r = 0;

  struct fem *fem = (struct fem *)sim->slv;

  struct vec w0;
  struct vec w1;
  struct vec w2;
  struct vec w3;

  if ((r = vec_new(&fem->prv.tmp, fem->prv.vec.n)))
    goto end;

  if ((r = vec_new(&w0, fem->prv.vec.n)))
    goto end;

  if ((r = vec_new(&w1, fem->prv.vec.n)))
    goto end;

  fem->slv.run.tb = 2;

  if (fem->ops.tdd > 2) {
    if ((r = vec_new(&w2, fem->prv.vec.n)))
      goto end;

    fem->slv.run.tb = 3;
  }

  if (fem->ops.tdd > 3) {
    if ((r = vec_new(&w3, fem->prv.vec.n)))
      goto end;

    fem->slv.run.tb = 4;
  }

  fem->slv.run.wgt[0] = &w0;
  fem->slv.run.wgt[1] = &w1;
  fem->slv.run.wgt[2] = &w2;
  fem->slv.run.wgt[3] = &w3;

  double beg = sim->ops.tdd.beg;
  double hop = sim->ops.tdd.hop;
  int num = sim->ops.tdd.num;

  fem->slv.run.ti = 0;
  fem->slv.run.tv = beg;

  if ((r = pbc_i1s_slv(sim)))
    goto end;

  for (int i = 1, t = hop; i <= num; i += 1, t += hop) {
    fem->slv.run.ti = i;
    fem->slv.run.tv = t;

    if ((r = sys_slv(sim)))
      goto end;

    if (fem->slv.itr.run)
      fem->slv.itr.run(fem->slv.itr.ctx, sim);

    if (sim->ops.exp.put)
      sim->ops.exp.put(sim);

    double *tmp = 0;

    switch (fem->slv.run.tb) {
    case 2:
      tmp = w1.dat;
      w1.dat = w0.dat;
      w0.dat = tmp;
      break;
    case 3:
      tmp = w2.dat;
      w2.dat = w1.dat;
      w1.dat = w0.dat;
      w0.dat = tmp;
      break;
    case 4:
      tmp = w3.dat;
      w3.dat = w2.dat;
      w2.dat = w1.dat;
      w1.dat = w0.dat;
      w0.dat = tmp;
      break;
    }
  }

end:
  vec_cls(&w0);
  vec_cls(&w1);
  vec_cls(&w2);
  vec_cls(&w3);

  return r;
}

static int pbc_i1s_slv(struct sim *sim) {
  struct fem *fem = (struct fem *)sim->slv;

  for (int i = 0; i < sim->msh->hxd.len; ++i) {
    struct hxd *hxd = &sim->msh->hxd.dat[i];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct cnd_ini *ini = &sim->cnd_ini.dat[obj->ini];

    if (ini->tgt.type == VAL_FUN) {
      struct sim_fun_ctx ctx = {
          .sim = sim,
          .vtx = -1,
          .qud = -1,
          .hxd = i,
      };

      for (int k = 0; k < 8; ++k) {
        ctx.vtx = hxd->vtx[k];
        fem->slv.run.wgt[1]->dat[hxd->vtx[k]] =
            ini->tgt.as.fun(&ctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
      }
    } else {
      for (int k = 0; k < 8; ++k)
        fem->slv.run.wgt[1]->dat[hxd->vtx[k]] = ini->tgt.as.num;
    }
  }

  return 0;
}

static int slv_fpi(struct sim *sim);
static int slv_new(struct sim *sim);

static int sys_slv(struct sim *sim) {
  assert(sim);

  struct fem *fem = (struct fem *)sim->slv;

  if (fem->ops.non.ops.fd) {
    switch (fem->ops.non.mod) {
    case NON_FPI:
      return slv_fpi(sim);
    case NON_NEW:
      return slv_new(sim);
    }
  }

  if (fem_lin_asm(sim))
    return -1;

  switch (fem->ops.iss.mod) {
  case ISS_BCG:
    if (iss_bcg_slv(&fem->prv.mtx, fem->slv.run.wgt[0], &fem->prv.vec,
                    &fem->ops.iss.ops.bcg))
      return -1;

    break;
  default:
    errno = ENOTSUP;
    return -1;
  }

  return 0;
}

static int slv_fpi(struct sim *sim) {
  int r = 0;

  struct fem *fem = (struct fem *)sim->slv;
  struct vec tmp;

  if ((r = vec_new(&tmp, fem->prv.vec.n)))
    goto end;

  int max = fem->ops.non.ops.max;
  double err = fem->ops.non.ops.err;
  double cur = 1;
  double tnm = 0;
  double bnm = 0;

  for (int i = 0; i < fem->prv.vec.n; ++i)
    fem->slv.run.wgt[0]->dat[i] = 1;

  for (int i = 0; i < max && cur > err; ++i) {
    if (i > 0)
      switch (fem->ops.iss.mod) {
      case ISS_BCG:
        if ((r = iss_bcg_slv(&fem->prv.mtx, fem->slv.run.wgt[0], &fem->prv.vec,
                             &fem->ops.iss.ops.bcg)))
          goto end;

        break;
      default:
        errno = ENOTSUP;
        r = -1;
        goto end;
      }

    smtx_rst(&fem->prv.mtx);
    // smtx_rst(&fem->prv.sig);
    // smtx_rst(&fem->prv.chi);
    vec_rst(&fem->prv.vec);

    if ((r = fem_lin_asm(sim)))
      goto end;

    mtx_vmlt(&fem->prv.mtx, fem->slv.run.wgt[0], &tmp);
    vec_cmb(&tmp, &fem->prv.vec, &tmp, -1);
    vec_nrm(&tmp, &tnm);
    vec_nrm(&fem->prv.vec, &bnm);

    cur = tnm / bnm;

    fem->ops.non.run.itr = i;
    fem->ops.non.run.err = cur;

    if (fem->ops.non.ops.itr.run)
      fem->ops.non.ops.itr.run(fem->ops.non.ops.itr.ctx, sim);
  }

end:
  vec_cls(&tmp);

  return r;
}

static int slv_new_asm(struct sim *sim);

static int slv_new(struct sim *sim) {
  (void)sim;
  int r = 0;

  struct fem *fem = (struct fem *)sim->slv;
  struct vec tmp;

  if ((r = vec_new(&tmp, fem->prv.vec.n)))
    goto end;

  int max = fem->ops.non.ops.max;
  double err = fem->ops.non.ops.err;
  double cur = 1;
  double tnm = 0;
  double bnm = 0;

  for (int i = 0; i < fem->prv.vec.n; ++i)
    fem->slv.run.wgt[0]->dat[i] = 1;

  for (int i = 0; i < max && cur > err; ++i) {
    if (i > 0) {
      if (slv_new_asm(sim))
        return -1;

      switch (fem->ops.iss.mod) {
      case ISS_BCG:
        if ((r = iss_bcg_slv(&fem->prv.mtx, fem->slv.run.wgt[0], &fem->prv.vec,
                             &fem->ops.iss.ops.bcg)))
          goto end;

        break;
      default:
        errno = ENOTSUP;
        r = -1;
        goto end;
      }
    }

    smtx_rst(&fem->prv.mtx);
    vec_rst(&fem->prv.vec);

    if ((r = fem_lin_asm(sim)))
      goto end;

    mtx_vmlt(&fem->prv.mtx, fem->slv.run.wgt[0], &tmp);
    vec_cmb(&tmp, &fem->prv.vec, &tmp, -1);
    vec_nrm(&tmp, &tnm);
    vec_nrm(&fem->prv.vec, &bnm);

    cur = tnm / bnm;

    fem->ops.non.run.itr = i;
    fem->ops.non.run.err = cur;

    if (fem->ops.non.ops.itr.run)
      fem->ops.non.ops.itr.run(fem->ops.non.ops.itr.ctx, sim);
  }

end:
  vec_cls(&tmp);

  return r;

  return 0;
}

double pdif(struct sim_fun_ctx *ctx, double (*fun)(void *, struct vec *),
            int var, double hop, struct vec *vtx) {

  ctx->sim->slv->run.wgt[0]->dat[var] += hop;
  double f1 = fun(ctx, vtx);

  ctx->sim->slv->run.wgt[0]->dat[var] -= 2 * hop;
  double f2 = fun(ctx, vtx);

  ctx->sim->slv->run.wgt[0]->dat[var] += hop;

  return (f1 - f2) / (2 * hop);
}

static int slv_new_asm(struct sim *sim) {
  static double gnx[2][2];
  static double gny[2][2];
  static double gnz[2][2];

  static double mnx[2][2][2];
  static double mny[2][2][2];
  static double mnz[2][2][2];

  struct fem *fem = (struct fem *)sim->slv;
  struct vec *vtx = sim->msh->vtx.dat;

  struct sim_fun_ctx ctx = {.sim = sim};

  for (int h = 0; h < sim->msh->hxd.len; ++h) {
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    int v0 = hxd->vtx[0];
    int v7 = hxd->vtx[7];

    double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
    double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
    double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) {
        gnx[i][j] = GN[i][j] / hx;
        gny[i][j] = GN[i][j] / hy;
        gnz[i][j] = GN[i][j] / hz;

        for (int k = 0; k < 2; ++k) {
          mnx[i][j][k] = MN[i][j][k] * hx;
          mny[i][j][k] = MN[i][j][k] * hy;
          mnz[i][j][k] = MN[i][j][k] * hz;
        }
      }

    for (int i = 0; i < 8; ++i) {
      int gi = hxd->vtx[i];

      int mui = MU[i];
      int nui = NU[i];
      int tti = TT[i];

      double bi = 0;

      for (int j = 0; j < 8; ++j) {
        int gj = hxd->vtx[j];

        int muj = MU[j];
        int nuj = NU[j];
        int ttj = TT[j];

        double mij = 0;
        double bij = 0;

        for (int k = 0; k < 8; ++k) {
          int gk = hxd->vtx[k];

          int muk = MU[k];
          int nuk = NU[k];
          int ttk = TT[k];

          mij += sim->slv->run.wgt[0]->dat[gk] *
                 pdif(&ctx, (double (*)(void *, struct vec *))mat->lam.as.fun,
                      j, HOP, &vtx[gj]) *
                 (gnx[muk][mui] * mny[nuj][nuk][nui] * mnz[ttj][ttk][tti] +
                  mnx[muj][muk][mui] * gny[nuk][nui] * mnz[ttj][ttk][tti] +
                  mnx[muj][muk][mui] * mny[nuj][nuk][nui] * gnz[ttk][tti]);

          bij += sim->slv->run.wgt[0]->dat[gk] *
                 pdif(&ctx, (double (*)(void *, struct vec *))mat->lam.as.fun,
                      k, HOP, &vtx[gk]) *
                 (gnx[muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * gny[nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * mny[nuk][nuj][nui] * gnz[ttj][tti]);
        }

        bi += bij * sim->slv->run.wgt[0]->dat[gj];

        if (smtx_add(&fem->prv.mtx, gi, gj, mij))
          return -1;
      }

      fem->prv.vec.dat[gi] += bi;
    }
  }

  return 0;
}
