#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/non/dif.h>
#include <numx/pde/sim.h>

#include "fem.h"

#define min(a, b) (((a) < (b)) ? (a) : (b))

static const double C = 10e10;

static const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
static const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
static const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

static const double G[2][2] = {
    {1.0,  -1.0},
    {-1.0, 1.0 }
};

static const double M[2][2] = {
    {2.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 6.0}
};

static const double GN[2][2] = {
    {1.0 / 2.0,  -1.0 / 2.0},
    {-1.0 / 2.0, 1.0 / 2.0 }
};

static const double MN[2][2][2] = {
    {{1.0 / 4.0, 1.0 / 12.0},  {1.0 / 12.0, 1.0 / 12.0}},
    {{1.0 / 12.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 4.0} }
};

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

static int ell_asm(struct sim *sim, struct fem_std_ctx *ctx);
static int pbc_asm(struct sim *sim, struct fem_std_ctx *ctx);
static int hyp_asm(struct sim *sim, struct fem_std_ctx *ctx);

int fem_std_lin_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    assert(sim);
    assert(ctx);

    switch (sim->eqn) {
        case SIM_ELL:
            return ell_asm(sim, ctx);
        case SIM_PBC:
            return pbc_asm(sim, ctx);
        case SIM_HYP:
            return hyp_asm(sim, ctx);
    }

    return 0;
}

struct asm_ops
{
    struct smtx *mlam;
    struct smtx *mgam;
    struct smtx *msig;
    struct smtx *mchi;
    struct smtx *mdir;
    struct smtx *mrob;

    struct vec *vsrc;
    struct vec *vdir;
    struct vec *vneu;
    struct vec *vrob;
};

static int assemble(struct sim *sim, struct asm_ops ops);

static int ell_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    mtx_rst(&ctx->mtx);
    vec_rst(&ctx->vec);

    assemble(sim, (struct asm_ops){
                      .mlam = &ctx->mtx,
                      .mgam = &ctx->mtx,
                      .msig = NULL,
                      .mchi = NULL,
                      .mdir = &ctx->mtx,
                      .mrob = &ctx->mtx,
                      .vsrc = &ctx->vec,
                      .vdir = &ctx->vec,
                      .vneu = &ctx->vec,
                      .vrob = &ctx->vec,
                  });

    return 0;
}

static int pbc_asm_i2s(struct sim *sim, struct fem_std_ctx *ctx);
static int pbc_asm_i3s(struct sim *sim, struct fem_std_ctx *ctx);
static int pbc_asm_i4s(struct sim *sim, struct fem_std_ctx *ctx);

static int pbc_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    int itr = sim->slv->run.ti;

    mtx_rst(&ctx->mtx);
    mtx_rst(&ctx->sig);
    vec_rst(&ctx->vec);

    assemble(sim, (struct asm_ops){
                      .mlam = &ctx->mtx,
                      .mgam = &ctx->mtx,
                      .msig = &ctx->sig,
                      .mchi = NULL,
                      .mdir = NULL,
                      .mrob = NULL,
                      .vsrc = NULL,
                      .vdir = NULL,
                      .vneu = NULL,
                      .vrob = NULL,
                  });

    static int (*f[3])(struct sim *, struct fem_std_ctx *) = {
        pbc_asm_i2s,
        pbc_asm_i3s,
        pbc_asm_i4s,
    };

    f[min(itr, (int)sim->slv->ops.tdd - 1) - 1](sim, ctx);

    return 0;
}

static int pbc_asm_i2s(struct sim *sim, struct fem_std_ctx *ctx)
{
    // int itr = sim->slv->run.ti;
    double hop = sim->ops.tdd.hop;

    mtx_cmb(&ctx->mtx, &ctx->sig, &ctx->mtx, 1.0 / hop);

    assemble(sim, (struct asm_ops){
                      .mlam = NULL,
                      .mgam = NULL,
                      .msig = NULL,
                      .mchi = NULL,
                      .mdir = &ctx->mtx,
                      .mrob = &ctx->mtx,
                      .vsrc = &ctx->vec,
                      .vdir = &ctx->vec,
                      .vneu = &ctx->vec,
                      .vrob = &ctx->vec,
                  });

    mtx_vmul(&ctx->sig, &ctx->w1, &ctx->tmp);
    vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, 1.0 / hop);

    return 0;
}

static int pbc_asm_i3s(struct sim *sim, struct fem_std_ctx *ctx)
{
    // int itr = sim->slv->run.ti;
    double hop = sim->ops.tdd.hop;

    mtx_cmb(&ctx->mtx, &ctx->sig, &ctx->mtx, 3.0 / (2 * hop));

    assemble(sim, (struct asm_ops){
                      .mlam = NULL,
                      .mgam = NULL,
                      .msig = NULL,
                      .mchi = NULL,
                      .mdir = &ctx->mtx,
                      .mrob = &ctx->mtx,
                      .vsrc = &ctx->vec,
                      .vdir = &ctx->vec,
                      .vneu = &ctx->vec,
                      .vrob = &ctx->vec,
                  });

    if (mtx_vmul(&ctx->sig, &ctx->w1, &ctx->tmp))
        return -1;

    if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, 2.0 / hop))
        return -1;

    if (mtx_vmul(&ctx->sig, &ctx->w2, &ctx->tmp))
        return -1;

    if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, -1.0 / (2 * hop)))
        return -1;

    return 0;
}

static int pbc_asm_i4s(struct sim *sim, struct fem_std_ctx *ctx)
{
    // int itr = sim->slv->run.ti;
    double hop = sim->ops.tdd.hop;

    mtx_cmb(&ctx->mtx, &ctx->sig, &ctx->mtx, 11.0 / (6 * hop));

    assemble(sim, (struct asm_ops){
                      .mlam = NULL,
                      .mgam = NULL,
                      .msig = NULL,
                      .mchi = NULL,
                      .mdir = &ctx->mtx,
                      .mrob = &ctx->mtx,
                      .vsrc = &ctx->vec,
                      .vdir = &ctx->vec,
                      .vneu = &ctx->vec,
                      .vrob = &ctx->vec,
                  });

    if (mtx_vmul(&ctx->sig, &ctx->w1, &ctx->tmp))
        return -1;

    if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, 3.0 / hop))
        return -1;

    if (mtx_vmul(&ctx->sig, &ctx->w2, &ctx->tmp))
        return -1;

    if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, -3.0 / (2 * hop)))
        return -1;

    if (mtx_vmul(&ctx->sig, &ctx->w3, &ctx->tmp))
        return -1;

    if (vec_cmb(&ctx->vec, &ctx->tmp, &ctx->vec, 1.0 / (3 * hop)))
        return -1;

    return 0;
}

static int hyp_asm(struct sim *, struct fem_std_ctx *)
{
    exit(-1);
}

static int asm_qud_dim(
    struct sim *sim, struct qud *qud, double *hxi, double *hzt)
{
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

static int assemble(struct sim *sim, struct asm_ops ops)
{
    struct vec *vtx = sim->msh->vtx.dat;

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    if (ops.mlam || ops.mgam || ops.msig || ops.mchi || ops.vsrc) {
        for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
            static double lam[8];
            static double gam[8];
            static double sig[8];
            static double chi[8];
            static double src[8];

            struct hxd *hxd = &sim->msh->hxd.dat[hi];
            struct obj *obj = &sim->obj.dat[hxd->pid];
            struct mat *mat = &sim->mat.dat[obj->mat];
            struct val *val = &sim->src.dat[obj->src];

            double lam_n = mat->lam.as.num;
            double gam_n = mat->gam.as.num;
            double sig_n = mat->sig.as.num;
            double chi_n = mat->chi.as.num;
            double src_n = val->as.num;

            bool lam_num = ops.mlam && mat->lam.type == VAL_NUM;
            bool gam_num = ops.mgam && mat->gam.type == VAL_NUM;
            bool sig_num = ops.msig && mat->sig.type == VAL_NUM;
            bool chi_num = ops.mchi && mat->chi.type == VAL_NUM;
            bool src_num = ops.vsrc && val->type == VAL_NUM;

            bool lam_fun = ops.mlam && mat->lam.type == VAL_FUN;
            bool gam_fun = ops.mgam && mat->gam.type == VAL_FUN;
            bool sig_fun = ops.msig && mat->sig.type == VAL_FUN;
            bool chi_fun = ops.mchi && mat->chi.type == VAL_FUN;
            bool src_fun = ops.vsrc && val->type == VAL_FUN;

            bool has_fun = lam_fun || gam_fun || sig_fun || chi_fun;

            fun_ctx.hxd = hi;

            if (has_fun || src_fun)
                for (int k = 0; k < 8; ++k) {
                    fun_ctx.vtx = hxd->vtx[k];

                    struct vec *v = &vtx[fun_ctx.vtx];

                    if (lam_fun)
                        lam[k] = mat->lam.as.fun(&fun_ctx, v);

                    if (gam_fun)
                        gam[k] = mat->gam.as.fun(&fun_ctx, v);

                    if (sig_fun)
                        sig[k] = mat->sig.as.fun(&fun_ctx, v);

                    if (chi_fun)
                        chi[k] = mat->chi.as.fun(&fun_ctx, v);

                    if (src_fun)
                        src[k] = val->as.fun(&fun_ctx, v);
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

                double bi = 0;

                if (src_num)
                    bi = src_n * hx * hy * hz / 8;

                for (int j = 0; j < 8; ++j) {
                    int gj = hxd->vtx[j];

                    int muj = MU[j];
                    int nuj = NU[j];
                    int ttj = TT[j];

                    double glp = gx[muj][mui] *
                                 my[nuj][nui] *
                                 mz[ttj][tti] +
                                 mx[muj][mui] *
                                 gy[nuj][nui] *
                                 mz[ttj][tti] +
                                 mx[muj][mui] *
                                 my[nuj][nui] *
                                 gz[ttj][tti];

                    double mlp = mx[muj][mui] * my[nuj][nui] * mz[ttj][tti];

                    if (src_fun)
                        bi += src[j] * mlp;

                    if (lam_num)
                        mtx_inc(ops.mlam, gi, gj, lam_n * glp);

                    if (gam_num)
                        mtx_inc(ops.mgam, gi, gj, gam_n * mlp);

                    if (sig_num)
                        mtx_inc(ops.msig, gi, gj, sig_n * mlp);

                    if (chi_num)
                        mtx_inc(ops.mchi, gi, gj, chi_n * mlp);

                    if (has_fun) {
                        double mij_lam = 0;
                        double mij_gam = 0;
                        double mij_sig = 0;
                        double mij_chi = 0;

                        for (int k = 0; k < 8; ++k) {
                            int muk = MU[k];
                            int nuk = NU[k];
                            int ttk = TT[k];

                            double gln = gnx[muj][mui] *
                                         mny[nuk][nuj][nui] *
                                         mnz[ttk][ttj][tti] +
                                         mnx[muk][muj][mui] *
                                         gny[nuj][nui] *
                                         mnz[ttk][ttj][tti] +
                                         mnx[muk][muj][mui] *
                                         mny[nuk][nuj][nui] *
                                         gnz[ttj][tti];

                            double mln = mnx[muk][muj][mui] *
                                         mny[nuk][nuj][nui] *
                                         mnz[ttk][ttj][tti];

                            if (lam_fun)
                                mij_lam += lam[k] * gln;

                            if (gam_fun)
                                mij_gam += gam[k] * mln;

                            if (sig_fun)
                                mij_sig += sig[k] * mln;

                            if (chi_fun)
                                mij_chi += chi[k] * mln;
                        }

                        if (lam_fun)
                            mtx_inc(ops.mlam, gi, gj, mij_lam);

                        if (gam_fun)
                            mtx_inc(ops.mgam, gi, gj, mij_gam);

                        if (sig_fun)
                            mtx_inc(ops.msig, gi, gj, mij_sig);

                        if (chi_fun)
                            mtx_inc(ops.mchi, gi, gj, mij_chi);
                    }
                }

                if (ops.vsrc)
                    ops.vsrc->dat[gi] += bi;
            }
        }
    }

    if (ops.mrob || ops.vneu || ops.vrob || ops.mdir || ops.vdir) {
        struct ilog dir;

        if (log_new(&dir))
            return -1;

        for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
            struct qud     *qud = &sim->msh->qud.dat[qi];
            struct bnd     *bnd = &sim->bnd.dat[qud->pid];
            struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

            fun_ctx.qud = qi;
            fun_ctx.hxd = qud->hxd;

            double hxi = 0;
            double hzt = 0;

            asm_qud_dim(sim, qud, &hxi, &hzt);

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) {
                    mx[i][j] = M[i][j] * hxi;
                    mz[i][j] = M[i][j] * hzt;
                }

            switch (cnd->type) {
                case CND_BND_DIR:
                    if (log_add(&dir, qi))
                        return -1;

                    break;
                case CND_BND_NEU: {
                    if (ops.vneu == NULL)
                        continue;

                    struct val   *tta_val = &cnd->pps.neu.tta;
                    static double tta[4];

                    double tta_n = tta_val->as.num;

                    if (tta_val->type == VAL_FUN)
                        for (int k = 0; k < 4; ++k) {
                            fun_ctx.vtx = qud->vtx[k];

                            tta[k] =
                                tta_val->as.fun(&fun_ctx, &vtx[fun_ctx.vtx]);
                        }

                    for (int i = 0; i < 4; ++i) {
                        int gi = qud->vtx[i];

                        int mui = MU[i];
                        int nui = NU[i];

                        double bi = 0;

                        if (tta_val->type == VAL_FUN) {
                            for (int k = 0; k < 4; ++k) {
                                int muk = MU[k];
                                int nuk = NU[k];

                                bi += tta[k] * mx[muk][mui] * mz[nuk][nui];
                            }
                        } else {
                            bi = tta_n * hxi * hzt / 4;
                        }

                        ops.vneu->dat[gi] += bi;
                    }

                    break;
                }
                case CND_BND_ROB: {
                    if (ops.vrob == NULL && ops.mrob == NULL)
                        continue;

                    struct val *bet_val = &cnd->pps.rob.bet;
                    struct val *ext_val = &cnd->pps.rob.ext;

                    double bet_n = cnd->pps.rob.bet.as.num;
                    double ext_n = cnd->pps.rob.ext.as.num;

                    bool bet_fun = bet_val->type == VAL_FUN;
                    bool ext_fun = ops.mrob && ext_val->type == VAL_FUN;

                    static double bet[4];
                    static double ext[4];

                    if (bet_fun || ext_fun)
                        for (int k = 0; k < 4; ++k) {
                            fun_ctx.vtx = qud->vtx[k];

                            if (bet_fun)
                                bet[k] = bet_val->as.fun(
                                    &fun_ctx, &vtx[fun_ctx.vtx]);

                            if (ext_fun)
                                ext[k] = ext_val->as.fun(
                                    &fun_ctx, &vtx[fun_ctx.vtx]);
                        }

                    switch (cnd->pps.rob.bet.type) {
                        case VAL_NUM:
                            if (ops.mrob)
                                for (int i = 0; i < 4; ++i) {
                                    int gi = qud->vtx[i];

                                    int mui = MU[i];
                                    int nui = NU[i];

                                    for (int j = 0; j < 4; ++j) {
                                        int gj = qud->vtx[j];

                                        int muj = MU[j];
                                        int nuj = NU[j];

                                        double mij = bet_n * (mx[muj][mui] *
                                                                 mz[nuj][nui]);

                                        mtx_inc(ops.mrob, gi, gj, mij);
                                    }
                                }

                            if (ops.vrob)
                                switch (ext_val->type) {
                                    case VAL_NUM:
                                        for (int i = 0; i < 4; ++i)
                                            ops.vrob->dat[qud->vtx[i]] +=
                                                bet_n * ext_n * hxi * hzt / 4;

                                        break;
                                    case VAL_FUN:
                                        for (int i = 0; i < 4; ++i) {
                                            int mui = MU[i];
                                            int nui = NU[i];

                                            double bi = 0;

                                            for (int k = 0; k < 4; ++k) {
                                                int muk = MU[k];
                                                int nuk = NU[k];

                                                bi += ext[k] *
                                                      mx[muk][mui] *
                                                      mz[nuk][nui];
                                            }

                                            ops.vrob->dat[qud->vtx[i]] +=
                                                bet_n * bi;
                                        }

                                        break;
                                    case VAL_HMC:
                                        return -1;
                                }

                            break;
                        case VAL_FUN:
                            if (ops.mrob)
                                for (int i = 0; i < 4; ++i) {
                                    int gi = qud->vtx[i];

                                    int mui = MU[i];
                                    int nui = NU[i];

                                    for (int j = 0; j < 4; ++j) {
                                        int gj = qud->vtx[j];

                                        int muj = MU[j];
                                        int nuj = NU[j];

                                        double mij = 0;

                                        for (int k = 0; k < 4; ++k) {
                                            int muk = MU[k];
                                            int nuk = NU[k];

                                            mij += bet[k] *
                                                   (mnx[muk][muj][mui] *
                                                       mnz[nuk][nuj][nui]);
                                        }

                                        mtx_inc(ops.mrob, gi, gj, mij);
                                    }
                                }

                            if (ops.vrob)
                                switch (cnd->pps.rob.ext.type) {
                                    case VAL_NUM:
                                        for (int i = 0; i < 4; ++i) {
                                            int mui = MU[i];
                                            int nui = NU[i];

                                            double bi = 0;

                                            for (int k = 0; k < 4; ++k) {
                                                int muk = MU[k];
                                                int nuk = NU[k];

                                                bi += bet[k] *
                                                      mx[muk][mui] *
                                                      mz[nuk][nui];
                                            }

                                            ops.vrob->dat[qud->vtx[i]] +=
                                                ext_n * bi;
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

                                                    bi += bet[k] *
                                                          ext[j] *
                                                          mnx[muk][muj][mui] *
                                                          mnz[nuk][nuj][nui];
                                                }
                                            }

                                            ops.vrob->dat[qud->vtx[i]] += bi;
                                        }

                                        break;
                                    case VAL_HMC:
                                        return -1;
                                }

                            break;
                        case VAL_HMC:
                            return -1;
                    }

                    break;
                }
            }
        }

        log_rst(&dir);

        if (ops.vdir || ops.mdir) {
            for (int qi = 0; !log_adv(&dir, &qi);) {
                struct qud     *qud = &sim->msh->qud.dat[qi];
                struct bnd     *bnd = &sim->bnd.dat[qud->pid];
                struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];
                struct val     *tgt_val = &cnd->pps.dir.tgt;

                fun_ctx.qud = qi;
                fun_ctx.hxd = qud->hxd;

                double tgt_n = tgt_val->as.num;

                for (int i = 0; i < 4; ++i) {
                    int gi = qud->vtx[i];

                    if (ops.mdir)
                        ops.mdir->dr[qud->vtx[i]] = C;

                    fun_ctx.vtx = gi;

                    if (ops.vdir == NULL)
                        continue;

                    if (tgt_val->type == VAL_FUN)
                        ops.vdir->dat[gi] =
                            C * tgt_val->as.fun(&fun_ctx, &vtx[gi]);
                    else
                        ops.vdir->dat[gi] = C * tgt_n;
                }
            }

            if (errno != ENOENT)
                return -1;

            errno = 0;
        }
    }

    return 0;
}

static void dif_twk_wgt(void *ctx, double hop, struct dif_ops *ops)
{
    ((struct sim_fun_ctx *)ctx)->sim->slv->run.wgt[0]->dat[ops->var] += hop;
}

int fem_std_lin_new(struct sim *sim, struct fem_std_ctx *ctx)
{
    struct vec *vtx = sim->msh->vtx.dat;

    struct sim_fun_ctx fun_ctx = {.sim = sim, .vtx = -1, .qud = -1, .hxd = -1};
    struct dif_ops     dif_ops = {.hop = 0, .twk = dif_twk_wgt};

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

                    dlam[k] = dif_tpm(&fun_ctx,
                        (double (*)(void *, struct vec *))mat->lam.as.fun,
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

                    mij += ctx->w0.dat[gk] *
                           dlam[j] *
                           (gnx[muk][mui] *
                               mny[nuj][nuk][nui] *
                               mnz[ttj][ttk][tti] +
                               mnx[muj][muk][mui] *
                               gny[nuk][nui] *
                               mnz[ttj][ttk][tti] +
                               mnx[muj][muk][mui] *
                               mny[nuj][nuk][nui] *
                               gnz[ttk][tti]);

                    bij += ctx->w0.dat[gk] *
                           dlam[k] *
                           (gnx[muj][mui] *
                               mny[nuk][nuj][nui] *
                               mnz[ttk][ttj][tti] +
                               mnx[muk][muj][mui] *
                               gny[nuj][nui] *
                               mnz[ttk][ttj][tti] +
                               mnx[muk][muj][mui] *
                               mny[nuk][nuj][nui] *
                               gnz[ttj][tti]);
                }

                bi += bij * ctx->w0.dat[gj];

                mtx_inc(&ctx->mtx, gi, gj, mij);
            }

            ctx->vec.dat[gi] += bi;
        }
    }

    return 0;
}
