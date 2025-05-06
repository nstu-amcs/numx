#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/non/dif.h>
#include <numx/pde/sim.h>

#include <prv/fem/lin.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))

static const double C = 10e10;

static const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
static const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
static const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

static const double G[2][2] = {{1.0, -1.0}, {-1.0, 1.0}};

static const double M[2][2] = {{2.0 / 6.0, 1.0 / 6.0}, {1.0 / 6.0, 2.0 / 6.0}};

static const double GN[2][2] = {{1.0 / 2.0, -1.0 / 2.0},
                                {-1.0 / 2.0, 1.0 / 2.0}};

static const double MN[2][2][2] = {
    {{1.0 / 4.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 12.0}},
    {{1.0 / 12.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 4.0}}};

static double gx[2][2];
static double gy[2][2];
static double gz[2][2];

static double mx[2][2];
static double my[2][2];
static double mz[2][2];

static double gnx[2][2];
static double gny[2][2];
static double gnz[2][2];

static double mnx[2][2][2];
static double mny[2][2][2];
static double mnz[2][2][2];

static int ell_asm(struct sim *sim, struct fem_ctx *ctx);
static int pbc_asm(struct sim *sim, struct fem_ctx *ctx);
static int hyp_asm(struct sim *sim, struct fem_ctx *ctx);

int fem_lin_asm(struct sim *sim, struct fem_ctx *ctx) {
  assert(sim);
  assert(ctx);

  switch (sim->mod) {
  case SIM_ELL:
    return ell_asm(sim, ctx);
  case SIM_PBC:
    return pbc_asm(sim, ctx);
  case SIM_HYP:
    return hyp_asm(sim, ctx);
  }

  return 0;
}

static int mtx_lam_asm(struct sim *sim, struct smtx *mtx, int h);
static int mtx_gam_asm(struct sim *sim, struct smtx *mtx, int h);
static int mtx_sig_asm(struct sim *sim, struct smtx *mtx, int h);
static int mtx_chi_asm(struct sim *sim, struct smtx *mtx, int h);
static int mtx_rob_asm(struct sim *sim, struct smtx *mtx, int q);
static int mtx_dir_asm(struct sim *sim, struct smtx *mtx, int q);

static int vec_src_asm(struct sim *sim, struct vec *vec, int h);
static int vec_dir_asm(struct sim *sim, struct vec *vec, int q);
static int vec_neu_asm(struct sim *sim, struct vec *vec, int q);
static int vec_rob_asm(struct sim *sim, struct vec *vec, int q);

static int ell_asm(struct sim *sim, struct fem_ctx *ctx) {
  mtx_rst(&ctx->mtx);
  vec_rst(&ctx->vec);

  for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
    if (mtx_lam_asm(sim, &ctx->mtx, hi))
      return -1;

    if (mtx_gam_asm(sim, &ctx->mtx, hi))
      return -1;

    if (vec_src_asm(sim, &ctx->vec, hi))
      return -1;
  }

  struct ilog dir;

  if (log_new(&dir))
    return -1;

  for (int qi = 0; qi < sim->msh->qud.len; ++qi)
    switch (
        sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type) {
    case CND_BND_DIR:
      if (log_add(&dir, qi))
        return -1;

      break;
    case CND_BND_NEU:
      if (vec_neu_asm(sim, &ctx->vec, qi))
        return -1;

      break;
    case CND_BND_ROB:
      if (mtx_rob_asm(sim, &ctx->mtx, qi))
        return -1;

      if (vec_rob_asm(sim, &ctx->vec, qi))
        return -1;

      break;
    }

  log_rst(&dir);

  for (int qi = 0; !log_adv(&dir, &qi);) {
    if (mtx_dir_asm(sim, &ctx->mtx, qi))
      return -1;

    if (vec_dir_asm(sim, &ctx->vec, qi))
      return -1;
  }

  if (errno != ENOENT)
    return -1;

  errno = 0;
  return 0;
}

static int pbc_asm_i2s(struct sim *sim, struct fem_ctx *ctx);
static int pbc_asm_i3s(struct sim *sim, struct fem_ctx *ctx);
static int pbc_asm_i4s(struct sim *sim, struct fem_ctx *ctx);

static int pbc_asm(struct sim *sim, struct fem_ctx *ctx) {
  struct fem *fem = (struct fem *)sim->slv;

  vec_rst(&ctx->vec);

  if (sim->slv->run.ti == 1) {
    for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
      mtx_lam_asm(sim, &ctx->mtx, hi);
      mtx_gam_asm(sim, &ctx->mtx, hi);
      mtx_sig_asm(sim, &ctx->sig, hi);
    }
  }

  static int (*f[3])(struct sim *, struct fem_ctx *) = {
      pbc_asm_i2s,
      pbc_asm_i3s,
      pbc_asm_i4s,
  };

  int itr = sim->slv->run.ti;

  if (f[min(itr, (int)fem->slv.ops.tdd - 1) - 1](sim, ctx))
    return -1;

  return 0;
}

static int pbc_asm_i2s(struct sim *sim, struct fem_ctx *ctx) {
  int itr = sim->slv->run.ti;
  int hop = sim->ops.tdd.hop;

  if (itr == 1 && mtx_cmb(&ctx->mtx, &ctx->sig, &ctx->mtx, 1.0 / hop))
    return -1;

  for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
    if (vec_src_asm(sim, &ctx->vec, hi))
      return -1;

  struct ilog dir;

  if (log_new(&dir))
    return -1;

  for (int qi = 0; qi < sim->msh->qud.len; ++qi)
    switch (
        sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type) {
    case CND_BND_DIR:
      if (log_add(&dir, qi))
        return -1;

      break;
    case CND_BND_NEU:
      if (vec_neu_asm(sim, &ctx->vec, qi))
        return -1;

      break;
    case CND_BND_ROB:
      if (itr == 1 && mtx_rob_asm(sim, &ctx->mtx, qi))
        return -1;

      if (vec_rob_asm(sim, &ctx->vec, qi))
        return -1;

      break;
    }

  log_rst(&dir);

  for (int qi = 0; !log_adv(&dir, &qi);) {
    if (itr == 1 && mtx_dir_asm(sim, &ctx->mtx, qi))
      return -1;

    if (vec_dir_asm(sim, &ctx->vec, qi))
      return -1;
  }

  log_cls(&dir);

  if (errno != ENOENT)
    return -1;

  errno = 0;

  if (mtx_vmul(&ctx->sig, &ctx->w1, &ctx->tmp))
    return -1;

  if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, 1.0 / hop))
    return -1;

  return 0;
}

static int pbc_asm_i3s(struct sim *, struct fem_ctx *) { return 0; }

static int pbc_asm_i4s(struct sim *, struct fem_ctx *) { return 0; }

static int hyp_asm(struct sim *, struct fem_ctx *) { return 0; }

static int asm_qud_dim(struct sim *sim, struct qud *qud, double *hxi,
                       double *hzt) {
  double dat[3];

  struct vec nrm = {
      .n = 3,
      .dat = dat,
  };

  if (msh_qud_nrm(sim->msh, qud, &nrm))
    return -1;

  struct vec *vtx = sim->msh->vtx.dat;

  int v0 = qud->vtx[0];
  int v3 = qud->vtx[3];

  if (nrm.dat[0] != 0) {
    *hxi = vtx[v3].dat[1] - vtx[v0].dat[1];
    *hzt = vtx[v3].dat[2] - vtx[v0].dat[2];
  } else if (nrm.dat[1] != 0) {
    *hxi = vtx[v3].dat[0] - vtx[v0].dat[0];
    *hzt = vtx[v3].dat[2] - vtx[v0].dat[2];
  } else {
    *hxi = vtx[v3].dat[0] - vtx[v0].dat[0];
    *hzt = vtx[v3].dat[1] - vtx[v0].dat[1];
  }

  return 0;
}

static int mtx_lam_asm(struct sim *sim, struct smtx *mtx, int h) {
  static double lam[8];

  struct vec *vtx = sim->msh->vtx.dat;
  struct hxd *hxd = &sim->msh->hxd.dat[h];
  struct obj *obj = &sim->obj.dat[hxd->pid];
  struct mat *mat = &sim->mat.dat[obj->mat];

  double lam_n = mat->lam.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = -1,
      .hxd = h,
  };

  if (mat->lam.type == VAL_FUN)
    for (int k = 0; k < 8; ++k) {
      fun_ctx.vtx = hxd->vtx[k];
      lam[k] = mat->lam.as.fun(&fun_ctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
    }

  int v0 = hxd->vtx[0];
  int v7 = hxd->vtx[7];

  double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
  double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
  double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      gx[i][j] = G[i][j] / hx;
      gy[i][j] = G[i][j] / hy;
      gz[i][j] = G[i][j] / hz;

      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;

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

    for (int j = 0; j < 8; ++j) {
      int gj = hxd->vtx[j];

      int muj = MU[j];
      int nuj = NU[j];
      int ttj = TT[j];

      double mij = 0;

      if (mat->lam.type == VAL_FUN) {
        for (int k = 0; k < 8; ++k) {
          int muk = MU[k];
          int nuk = NU[k];
          int ttk = TT[k];

          mij += lam[k] *
                 (gnx[muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * gny[nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * mny[nuk][nuj][nui] * gnz[ttj][tti]);
        }
      } else {
        mij += lam_n * (gx[muj][mui] * my[nuj][nui] * mz[ttj][tti] +
                        mx[muj][mui] * gy[nuj][nui] * mz[ttj][tti] +
                        mx[muj][mui] * my[nuj][nui] * gz[ttj][tti]);
      }

      mtx_inc(mtx, gi, gj, mij);
    }
  }

  return 0;
}

static int mtx_gam_asm(struct sim *sim, struct smtx *mtx, int h) {
  static double gam[8];

  struct vec *vtx = sim->msh->vtx.dat;
  struct hxd *hxd = &sim->msh->hxd.dat[h];
  struct obj *obj = &sim->obj.dat[hxd->pid];
  struct mat *mat = &sim->mat.dat[obj->mat];

  double gam_n = mat->gam.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = -1,
      .hxd = h,
  };

  if (mat->gam.type == VAL_FUN)
    for (int k = 0; k < 8; ++k) {
      fun_ctx.vtx = hxd->vtx[k];
      gam[k] = mat->gam.as.fun(&fun_ctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
    }

  int v0 = hxd->vtx[0];
  int v7 = hxd->vtx[7];

  double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
  double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
  double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;

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

    for (int j = 0; j < 8; ++j) {
      int gj = hxd->vtx[j];

      int muj = MU[j];
      int nuj = NU[j];
      int ttj = TT[j];

      double mij = 0;

      if (mat->gam.type == VAL_FUN) {
        for (int k = 0; k < 8; ++k) {
          int muk = MU[k];
          int nuk = NU[k];
          int ttk = TT[k];

          mij += gam[k] *
                 (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
        }
      } else {
        mij += gam_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
      }

      mtx_inc(mtx, gi, gj, mij);
    }
  }

  return 0;
}

static int mtx_sig_asm(struct sim *sim, struct smtx *mtx, int h) {
  static double sig[8];

  struct vec *vtx = sim->msh->vtx.dat;
  struct hxd *hxd = &sim->msh->hxd.dat[h];
  struct obj *obj = &sim->obj.dat[hxd->pid];
  struct mat *mat = &sim->mat.dat[obj->mat];

  double sig_n = mat->sig.as.num;

  if (mat->sig.type == VAL_FUN)
    for (int k = 0; k < 8; ++k)
      sig[k] = mat->sig.as.fun(&((struct sim_fun_ctx){
                                   .sim = sim,
                                   .vtx = hxd->vtx[k],
                                   .qud = -1,
                                   .hxd = h,
                               }),
                               &sim->msh->vtx.dat[hxd->vtx[k]]);

  int v0 = hxd->vtx[0];
  int v7 = hxd->vtx[7];

  double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
  double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
  double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;

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

    for (int j = 0; j < 8; ++j) {
      int gj = hxd->vtx[j];

      int muj = MU[j];
      int nuj = NU[j];
      int ttj = TT[j];

      double mij = 0;

      if (mat->sig.type == VAL_FUN) {
        for (int k = 0; k < 8; ++k) {
          int muk = MU[k];
          int nuk = NU[k];
          int ttk = TT[k];

          mij += sig[k] *
                 (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
        }
      } else {
        mij += sig_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
      }

      if (mtx_inc(mtx, gi, gj, mij))
        return -1;
    }
  }

  return 0;
}

[[maybe_unused]]
static int mtx_chi_asm(struct sim *sim, struct smtx *mtx, int h) {
  static double chi[8];

  struct vec *vtx = sim->msh->vtx.dat;
  struct hxd *hxd = &sim->msh->hxd.dat[h];
  struct obj *obj = &sim->obj.dat[hxd->pid];
  struct mat *mat = &sim->mat.dat[obj->mat];

  double chi_n = mat->chi.as.num;

  if (mat->chi.type == VAL_FUN)
    for (int k = 0; k < 8; ++k)
      chi[k] = mat->chi.as.fun(&((struct sim_fun_ctx){
                                   .sim = sim,
                                   .vtx = hxd->vtx[k],
                                   .qud = -1,
                                   .hxd = h,
                               }),
                               &sim->msh->vtx.dat[hxd->vtx[k]]);

  int v0 = hxd->vtx[0];
  int v7 = hxd->vtx[7];

  double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
  double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
  double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;

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

    for (int j = 0; j < 8; ++j) {
      int gj = hxd->vtx[j];

      int muj = MU[j];
      int nuj = NU[j];
      int ttj = TT[j];

      double mij = 0;

      if (mat->chi.type == VAL_FUN) {
        for (int k = 0; k < 8; ++k) {
          int muk = MU[k];
          int nuk = NU[k];
          int ttk = TT[k];

          mij += chi[k] *
                 (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
        }
      } else {
        mij += chi_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
      }

      if (mtx_inc(mtx, gi, gj, mij))
        return -1;
    }
  }

  return 0;
}

static int mtx_rob_asm(struct sim *sim, struct smtx *mtx, int q) {
  static double bet[4];

  struct qud *qud = &sim->msh->qud.dat[q];
  struct bnd *bnd = &sim->bnd.dat[qud->pid];
  struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

  double bet_n = cnd->pps.rob.bet.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = q,
      .hxd = qud->hxd,
  };

  if (cnd->pps.rob.bet.type == VAL_FUN)
    for (int k = 0; k < 4; ++k) {
      fun_ctx.vtx = qud->vtx[k];
      bet[k] =
          cnd->pps.rob.bet.as.fun(&fun_ctx, &sim->msh->vtx.dat[qud->vtx[k]]);
    }

  double hxi = 0;
  double hzt = 0;

  if (asm_qud_dim(sim, qud, &hxi, &hzt))
    return -1;

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hxi;
      mz[i][j] = M[i][j] * hzt;

      for (int k = 0; k < 2; ++k) {
        mnx[i][j][k] = MN[i][j][k] * hxi;
        mnz[i][j][k] = MN[i][j][k] * hzt;
      }
    }

  for (int i = 0; i < 4; ++i) {
    int gi = qud->vtx[i];

    int mui = MU[i];
    int nui = NU[i];

    for (int j = 0; j < 4; ++j) {
      int gj = qud->vtx[j];

      int muj = MU[j];
      int nuj = NU[j];

      double mij = 0;

      switch (cnd->pps.rob.bet.type) {
      case VAL_NUM:
        mij = bet_n * (mx[muj][mui] * mz[nuj][nui]);
        break;
      case VAL_FUN:
        for (int k = 0; k < 4; ++k) {
          int muk = MU[k];
          int nuk = NU[k];

          mij += bet[k] * (mnx[muk][muj][mui] * mnz[nuk][nuj][nui]);
        }

        break;
      }

      mtx_inc(mtx, gi, gj, mij);
    }
  }

  return 0;
}

static int mtx_dir_asm(struct sim *sim, struct smtx *mtx, int q) {
  struct qud *qud = &sim->msh->qud.dat[q];

  for (int i = 0; i < 4; ++i)
    mtx->dr[qud->vtx[i]] = C;

  return 0;
}

static int vec_src_asm(struct sim *sim, struct vec *vec, int h) {
  static double src[8];

  struct vec *vtx = sim->msh->vtx.dat;
  struct hxd *hxd = &sim->msh->hxd.dat[h];
  struct obj *obj = &sim->obj.dat[hxd->pid];
  struct val *val = &sim->src.dat[obj->src];

  double src_n = val->as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = -1,
      .hxd = h,
  };

  if (val->type == VAL_FUN)
    for (int k = 0; k < 8; ++k) {
      fun_ctx.vtx = hxd->vtx[k];
      src[k] = val->as.fun(&fun_ctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
    }

  int v0 = hxd->vtx[0];
  int v7 = hxd->vtx[7];

  double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
  double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
  double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hx;
      my[i][j] = M[i][j] * hy;
      mz[i][j] = M[i][j] * hz;
    }

  for (int i = 0; i < 8; ++i) {
    int gi = hxd->vtx[i];

    int mui = MU[i];
    int nui = NU[i];
    int tti = TT[i];

    double bi = 0;

    if (val->type == VAL_FUN) {
      for (int k = 0; k < 8; ++k) {
        int muk = MU[k];
        int nuk = NU[k];
        int ttk = TT[k];

        bi += src[k] * (mx[muk][mui] * my[nuk][nui] * mz[ttk][tti]);
      }
    } else {
      bi = src_n * hx * hy * hz / 8;
    }

    vec->dat[gi] += bi;
  }

  return 0;
}

static int vec_dir_asm(struct sim *sim, struct vec *vec, int q) {
  struct qud *qud = &sim->msh->qud.dat[q];
  struct bnd *bnd = &sim->bnd.dat[qud->pid];
  struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

  double tgt_n = cnd->pps.dir.tgt.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = q,
      .hxd = qud->hxd,
  };

  for (int i = 0; i < 4; ++i) {
    int gi = qud->vtx[i];

    fun_ctx.vtx = gi;

    if (cnd->pps.dir.tgt.type == VAL_FUN)
      vec->dat[gi] =
          C * cnd->pps.dir.tgt.as.fun(&fun_ctx, &sim->msh->vtx.dat[gi]);
    else
      vec->dat[gi] = C * tgt_n;
  }

  return 0;
}

static int vec_neu_asm(struct sim *sim, struct vec *vec, int q) {
  static double tta[4];

  struct qud *qud = &sim->msh->qud.dat[q];
  struct bnd *bnd = &sim->bnd.dat[qud->pid];
  struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

  double tta_n = cnd->pps.neu.tta.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = q,
      .hxd = qud->hxd,
  };

  if (cnd->pps.neu.tta.type == VAL_FUN)
    for (int k = 0; k < 4; ++k) {
      fun_ctx.vtx = qud->vtx[k];
      tta[k] =
          cnd->pps.neu.tta.as.fun(&fun_ctx, &sim->msh->vtx.dat[fun_ctx.vtx]);
    }

  double hxi = 0;
  double hzt = 0;

  asm_qud_dim(sim, qud, &hxi, &hzt);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hxi;
      mz[i][j] = M[i][j] * hzt;
    }

  for (int i = 0; i < 4; ++i) {
    int gi = qud->vtx[i];

    int mui = MU[i];
    int nui = NU[i];

    double bi = 0;

    if (cnd->pps.neu.tta.type == VAL_FUN) {
      for (int k = 0; k < 4; ++k) {
        int muk = MU[k];
        int nuk = NU[k];

        bi += tta[k] * mx[muk][mui] * mz[nuk][nui];
      }
    } else {
      bi = tta_n * hxi * hzt / 4;
    }

    vec->dat[gi] += bi;
  }

  return 0;
}

static int vec_rob_asm(struct sim *sim, struct vec *vec, int q) {
  static double bet[4];
  static double ext[4];

  struct qud *qud = &sim->msh->qud.dat[q];
  struct bnd *bnd = &sim->bnd.dat[qud->pid];
  struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

  double bet_n = cnd->pps.rob.bet.as.num;
  double ext_n = cnd->pps.rob.ext.as.num;

  struct sim_fun_ctx fun_ctx = {
      .sim = sim,
      .vtx = -1,
      .qud = q,
      .hxd = qud->hxd,
  };

  if (cnd->pps.rob.bet.type == VAL_FUN)
    for (int k = 0; k < 4; ++k) {
      fun_ctx.vtx = qud->vtx[k];
      bet[k] =
          cnd->pps.rob.bet.as.fun(&fun_ctx, &sim->msh->vtx.dat[fun_ctx.vtx]);
    }

  if (cnd->pps.rob.ext.type == VAL_FUN)
    for (int k = 0; k < 4; ++k) {
      fun_ctx.vtx = qud->vtx[k];
      ext[k] =
          cnd->pps.rob.ext.as.fun(&fun_ctx, &sim->msh->vtx.dat[fun_ctx.vtx]);
    }

  double hxi = 0;
  double hzt = 0;

  if (asm_qud_dim(sim, qud, &hxi, &hzt))
    return -1;

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      mx[i][j] = M[i][j] * hxi;
      mz[i][j] = M[i][j] * hzt;

      for (int k = 0; k < 2; ++k) {
        mnx[i][j][k] = MN[i][j][k] * hxi;
        mnz[i][j][k] = MN[i][j][k] * hzt;
      }
    }

  switch (cnd->pps.rob.bet.type) {
  case VAL_NUM:
    switch (cnd->pps.rob.ext.type) {
    case VAL_NUM:
      for (int i = 0; i < 4; ++i)
        vec->dat[qud->vtx[i]] += bet_n * ext_n * hxi * hzt / 4;

      break;
    case VAL_FUN:
      for (int i = 0; i < 4; ++i) {
        int mui = MU[i];
        int nui = NU[i];

        double bi = 0;

        for (int k = 0; k < 4; ++k) {
          int muk = MU[k];
          int nuk = NU[k];

          bi += ext[k] * mx[muk][mui] * mz[nuk][nui];
        }

        vec->dat[qud->vtx[i]] += bet_n * bi;
      }

      break;
    }

    break;
  case VAL_FUN:
    switch (cnd->pps.rob.ext.type) {
    case VAL_NUM:
      for (int i = 0; i < 4; ++i) {
        int mui = MU[i];
        int nui = NU[i];

        double bi = 0;

        for (int k = 0; k < 4; ++k) {
          int muk = MU[k];
          int nuk = NU[k];

          bi += bet[k] * mx[muk][mui] * mz[nuk][nui];
        }

        vec->dat[qud->vtx[i]] += ext_n * bi;
      }

      break;
    case VAL_FUN:
      for (int i = 0; i < 4; ++i) {
        int mui = MU[i];
        int nui = NU[i];

        double bi = 0;

        for (int k = 0; k < 4; ++k) {
          int muk = MU[k];
          int nuk = NU[k];

          for (int j = 0; j < 4; ++j) {
            int muj = MU[j];
            int nuj = NU[j];

            bi += bet[k] * ext[j] * mnx[muk][muj][mui] * mnz[nuk][nuj][nui];
          }
        }

        vec->dat[qud->vtx[i]] += bi;
      }

      break;
    }

    break;
  }

  return 0;
}

static void dif_twk_wgt(void *ctx, double hop, struct dif_ops *ops) {
  ((struct sim_fun_ctx *)ctx)->sim->slv->run.wgt[0]->dat[ops->var] += hop;
}

int fem_lin_new(struct sim *sim, struct fem_ctx *ctx) {
  struct vec *vtx = sim->msh->vtx.dat;

  struct sim_fun_ctx fun_ctx = {.sim = sim, .vtx = -1, .qud = -1, .hxd = -1};
  struct dif_ops dif_ops = {.hop = 0, .twk = dif_twk_wgt};

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

    double dlam[8];

    switch (sim->slv->ops.non.ops.dif) {
    case DIF_GIV:
      for (int k = 0; k < 8; ++k) {
        fun_ctx.vtx = hxd->vtx[k];
        fun_ctx.hxd = h;
        dlam[k] = mat->lam.ops.dif(&fun_ctx, &vtx[fun_ctx.vtx]);
      }

      break;
    case DIF_NUM:
      for (int k = 0; k < 8; ++k) {
        int gk = hxd->vtx[k];

        dif_ops.var = gk;
        dif_ops.vtx = &vtx[gk];
        fun_ctx.vtx = gk;
        fun_ctx.hxd = h;

        dlam[k] =
            dif_tpm(&fun_ctx, (double (*)(void *, struct vec *))mat->lam.as.fun,
                    &dif_ops);
      }

      break;
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

          mij += ctx->w0.dat[gk] * dlam[j] *
                 (gnx[muk][mui] * mny[nuj][nuk][nui] * mnz[ttj][ttk][tti] +
                  mnx[muj][muk][mui] * gny[nuk][nui] * mnz[ttj][ttk][tti] +
                  mnx[muj][muk][mui] * mny[nuj][nuk][nui] * gnz[ttk][tti]);

          bij += ctx->w0.dat[gk] * dlam[k] *
                 (gnx[muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * gny[nuj][nui] * mnz[ttk][ttj][tti] +
                  mnx[muk][muj][mui] * mny[nuk][nuj][nui] * gnz[ttj][tti]);
        }

        bi += bij * ctx->w0.dat[gj];

        mtx_inc(&ctx->mtx, gi, gj, mij);
      }

      ctx->vec.dat[gi] += bi;
    }
  }

  return 0;
}
