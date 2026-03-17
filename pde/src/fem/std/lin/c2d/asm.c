#include <numx/com/log.h>
#include <stdio.h>

#include "../lin.h"

static double const DPH = 0.5;

/**
 * @brief Approximate solution at an arbitrary point.
 */
double fem_std_lin_c2d_apx(void *ctx, struct vec *vtx)
{
    assert(ctx);
    assert(vtx);

    struct apx_fun_ctx *apx_ctx = (struct apx_fun_ctx *)ctx;
    struct v2d         *v = apx_ctx->sim->msh->vtx.v2d.dat;

    if (apx_ctx->qud == -1) {
        apx_ctx->qud = umsh_vtx_qud_lup(apx_ctx->sim->msh, vtx);
    }

    if (apx_ctx->qud == -1) {
        return DBL_MAX;
    }

    if (apx_ctx->wgt == NULL) {
        apx_ctx->wgt = apx_ctx->sim->slv->run.wgt[0];
    }

    struct qud *q = &apx_ctx->sim->msh->qud.dat[apx_ctx->qud];

    double *w = apx_ctx->wgt->dat;

    int v0 = q->vtx[0];
    int v3 = q->vtx[3];

    double x1 = v[v0].dat[0];
    double x2 = v[v3].dat[0];
    double y1 = v[v0].dat[1];
    double y2 = v[v3].dat[1];

    double hm = (x2 - x1) * (y2 - y1);

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double r = 0;

    r += w[q->vtx[0]] * (x2 - x) * (y2 - y) / hm;
    r += w[q->vtx[1]] * (x - x1) * (y2 - y) / hm;
    r += w[q->vtx[2]] * (x2 - x) * (y - y1) / hm;
    r += w[q->vtx[3]] * (x - x1) * (y - y1) / hm;

    return r;
}

/**
 * @brief Approximate solution derivative at an arbitrary point.
 */
double fem_std_lin_c2d_dif(void *ctx, struct vec *vtx)
{
    assert(ctx);

    struct apx_fun_ctx *apx_ctx = (struct apx_fun_ctx *)ctx;
    struct v2d         *v = apx_ctx->sim->msh->vtx.v2d.dat;

    if (apx_ctx->qud == -1) {
        apx_ctx->qud = umsh_vtx_qud_lup(apx_ctx->sim->msh, vtx);
    }

    if (apx_ctx->qud == -1) {
        return DBL_MAX;
    }

    if (apx_ctx->wgt == NULL) {
        apx_ctx->wgt = apx_ctx->sim->slv->run.wgt[0];
    }

    struct qud *q = &apx_ctx->sim->msh->qud.dat[apx_ctx->qud];

    double *w = apx_ctx->wgt->dat;

    int v0 = q->vtx[0];
    int v3 = q->vtx[3];

    double x1 = v[v0].dat[0];
    double x2 = v[v3].dat[0];
    double y1 = v[v0].dat[1];
    double y2 = v[v3].dat[1];

    double hm = (x2 - x1) * (y2 - y1);

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double r = 0;

    switch (apx_ctx->var) {
        case 0: // x
            r += w[q->vtx[0]] * (-1) * (y2 - y) / hm;
            r += w[q->vtx[1]] * (1) * (y2 - y) / hm;
            r += w[q->vtx[2]] * (-1) * (y - y1) / hm;
            r += w[q->vtx[3]] * (1) * (y - y1) / hm;
            break;
        case 1: // y
            r += w[q->vtx[0]] * (x2 - x) * (-1) / hm;
            r += w[q->vtx[1]] * (x - x1) * (-1) / hm;
            r += w[q->vtx[2]] * (x2 - x) * (1) / hm;
            r += w[q->vtx[3]] * (x - x1) * (1) / hm;
            break;
        default:
            return DBL_MAX;
    }

    return r;
}

static int ell_asm(struct sim *sim, struct fem_std_ctx *ctx);
static int pbc_asm(struct sim *sim, struct fem_std_ctx *ctx);
static int hyp_asm(struct sim *sim, struct fem_std_ctx *ctx);

/**
 * @brief Assemble linear system for cortesian 2D mesh.
 */
int fem_std_lin_c2d_asm(struct sim *sim, struct fem_std_ctx *ctx)
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

    return -ENOTSUP;
}

struct asm_ops
{
    struct smtx *mlam; // Stiffness matrix.
    struct smtx *mgam; // Mass matrix.
    struct smtx *msig; // Parabolic matrix part.
    struct smtx *mchi; // Hyperbolic matrix part.
    struct smtx *mdir; // Dirichlet conditions matrix.
    struct smtx *mrob; // Robin conditions matrix.

    struct vec *vsrc; // Source vector.
    struct vec *vdir; // Dirichlet conditions vector.
    struct vec *vneu; // Neumann conditions vector.
    struct vec *vrob; // Robin conditions vector.
};

/**
 * @brief Assemble linear system using given options.
 */
static int assemble(struct sim *sim, struct asm_ops ops);

/**
 * @brief Assemble linear system for elliptic equation.
 */
static int ell_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    mtx_rst(&ctx->mtx);
    vec_rst(&ctx->vec);

    int r = assemble(sim, (struct asm_ops){
                              .mlam = &ctx->mtx,
                              .mgam = &ctx->mtx,
                              .msig = NULL, // not required
                              .mchi = NULL, // not required
                              .mdir = &ctx->mtx,
                              .mrob = &ctx->mtx,
                              .vsrc = &ctx->vec,
                              .vdir = &ctx->vec,
                              .vneu = &ctx->vec,
                              .vrob = &ctx->vec,
                          });

    if (!r) {
        printf("[fem][std][lin][c2d][asm] ok\n");
    }

    return r;
}

__attribute__((unused));
static int pbc_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    (void)sim;
    (void)ctx;

    return -ENOTSUP;
}

__attribute__((unused));
static int hyp_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    (void)sim;
    (void)ctx;

    return -ENOTSUP;
}

static const int MU[4] = {0, 1, 0, 1};
static const int NU[4] = {0, 0, 1, 1};

static const double M[2][2] = {
    {1.0 / 3.0, 1.0 / 6.0},
    {1.0 / 6.0, 1.0 / 3.0}
};

static const double GN[2][2][2] = {
    {{1.0 / 2.0, -1.0 / 2.0}, {-1.0 / 2.0, 1.0 / 2.0}},
    {{1.0 / 2.0, -1.0 / 2.0}, {-1.0 / 2.0, 1.0 / 2.0}},
};

static const double MN[2][2][2] = {
    {{1.0 / 4.0, 1.0 / 12.0},  {1.0 / 12.0, 1.0 / 12.0}},
    {{1.0 / 12.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 4.0} },
};

static double mx[2][2];
static double my[2][2];

static double gnx[2][2][2];
static double gny[2][2][2];
static double mnx[2][2][2];
static double mny[2][2][2];

static double F1(double a, double b, double v)
{
    return (b - v) / (b - a);
}

static double F2(double a, double b, double v)
{
    return (v - a) / (b - a);
}

static double (*F[2])(double a, double b, double v) = {F1, F2};

/**
 * @brief Get segment dimensions (size).
 */
static int asm_seg_dim(struct sim *sim, struct seg *seg, double *hxi)
{
    double     dat[2];
    struct vec nrm = {.n = 2, .dat = dat};

    if (umsh_seg_nrm(sim->msh, seg, &nrm)) {
        return -1;
    }

    struct v2d *vtx = sim->msh->vtx.v2d.dat;

    int v0 = seg->vtx[0];
    int v1 = seg->vtx[1];

    if (nrm.dat[0] != 0) { // x-norm -> x-fixed
        *hxi = vtx[v1].dat[1] - vtx[v0].dat[1];
        return 1;
    } else { // y-norm -> y-fixed
        *hxi = vtx[v1].dat[0] - vtx[v0].dat[0];
        return 2;
    }

    return -1;
}

static int asm_qud(struct sim *sim, struct asm_ops ops);
static int asm_seg(struct sim *sim, struct asm_ops ops);

static int assemble(struct sim *sim, struct asm_ops ops)
{

    int r = 0;

    if (ops.mlam || ops.mgam || ops.msig || ops.mchi || ops.vsrc) {
        if ((r = asm_qud(sim, ops))) {
            return r;
        }
    }

    if (ops.mrob || ops.vneu || ops.vrob || ops.mdir || ops.vdir) {
        if ((r = asm_seg(sim, ops))) {
            return r;
        }
    }

    return 0;
}

static int asm_qud(struct sim *sim, struct asm_ops ops)
{
    struct v2d *vtx = sim->msh->vtx.v2d.dat;

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        static double lam[4];
        static double gam[4];
        static double sig[4];
        static double chi[4];
        static double src[4];

        struct qud *qud = &sim->msh->qud.dat[qi];
        struct obj *obj = &sim->obj.dat[qud->pid];
        struct mat *mat = &sim->mat.dat[obj->mat];
        struct val *val = &sim->src.dat[obj->src];

        bool lam_fun = ops.mlam && mat->lam.type == VAL_FUN;
        bool gam_fun = ops.mgam && mat->gam.type == VAL_FUN;
        bool sig_fun = ops.msig && mat->sig.type == VAL_FUN;
        bool chi_fun = ops.mchi && mat->chi.type == VAL_FUN;
        bool src_fun = ops.vsrc && val->type == VAL_FUN;

        fun_ctx.qud = qi;

        // decompose parameters

        for (int k = 0; k < 4; ++k) {
            fun_ctx.vtx = qud->vtx[k];

            struct vec vw = {.dat = vtx[fun_ctx.vtx].dat, .n = 2};

            if (lam_fun) {
                fun_ctx.ctx = mat->lam.as.fun.ctx;
                lam[k] = mat->lam.as.fun.run(&fun_ctx, &vw);
            } else {
                lam[k] = mat->lam.as.num;
            }

            if (gam_fun) {
                fun_ctx.ctx = mat->gam.as.fun.ctx;
                gam[k] = mat->gam.as.fun.run(&fun_ctx, &vw);
            } else {
                gam[k] = mat->gam.as.num;
            }

            if (sig_fun) {
                fun_ctx.ctx = mat->sig.as.fun.ctx;
                sig[k] = mat->sig.as.fun.run(&fun_ctx, &vw);
            } else {
                sig[k] = mat->sig.as.num;
            }

            if (chi_fun) {
                fun_ctx.ctx = mat->chi.as.fun.ctx;
                chi[k] = mat->chi.as.fun.run(&fun_ctx, &vw);
            } else {
                chi[k] = mat->chi.as.num;
            }

            if (src_fun) {
                fun_ctx.ctx = val->as.fun.ctx;
                src[k] = val->as.fun.run(&fun_ctx, &vw);
            } else {
                src[k] = val->as.num;
            }
        }

        int v0 = qud->vtx[0];
        int v3 = qud->vtx[3];

        double hx = vtx[v3].dat[0] - vtx[v0].dat[0];
        double hy = vtx[v3].dat[1] - vtx[v0].dat[1];

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                mx[i][j] = M[i][j] * hx;
                my[i][j] = M[i][j] * hy;

                for (int k = 0; k < 2; ++k) {
                    gnx[i][j][k] = GN[i][j][k] / hx;
                    gny[i][j][k] = GN[i][j][k] / hy;
                    mnx[i][j][k] = MN[i][j][k] * hx;
                    mny[i][j][k] = MN[i][j][k] * hy;
                }
            }
        }

        for (int i = 0; i < 4; ++i) {
            int gi = qud->vtx[i]; // global i-index

            int mui = MU[i];
            int nui = NU[i];

            double bi = 0;

            for (int j = 0; j < 4; ++j) {
                int gj = qud->vtx[j]; // global j-index

                int muj = MU[j];
                int nuj = NU[j];

                bi += src[j] * mx[muj][mui] * my[nuj][nui];

                double mij_lam = 0;
                double mij_gam = 0;
                double mij_sig = 0;
                double mij_chi = 0;

                for (int k = 0; k < 4; ++k) {
                    int muk = MU[k];
                    int nuk = NU[k];

                    double gln = gnx[muk][muj][mui] *
                                 mny[nuk][nuj][nui] +
                                 mnx[muk][muj][mui] *
                                 gny[nuk][nuj][nui];
                    double mln = mnx[muk][muj][mui] * mny[nuk][nuj][nui];

                    mij_lam += lam[k] * gln;
                    mij_gam += gam[k] * mln;
                    mij_sig += sig[k] * mln;
                    mij_chi += chi[k] * mln;
                }

                if (ops.mlam) {
                    mtx_inc(ops.mlam, gi, gj, mij_lam);
                }

                if (ops.mgam) {
                    mtx_inc(ops.mgam, gi, gj, mij_gam);
                }

                if (ops.msig) {
                    mtx_inc(ops.msig, gi, gj, mij_sig);
                }

                if (ops.mchi) {
                    mtx_inc(ops.mchi, gi, gj, mij_chi);
                }
            }

            if (ops.vsrc) {
                ops.vsrc->dat[gi] += bi;
            }
        }
    }
    return 0;
}

static int asm_seg(struct sim *sim, struct asm_ops ops)
{
    struct v2d *vtx = sim->msh->vtx.v2d.dat;

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    for (int si = 0; si < sim->msh->seg.len; ++si) {
        struct seg     *seg = &sim->msh->seg.dat[si];
        struct qud     *qud = &sim->msh->qud.dat[seg->qud];
        struct bnd     *bnd = &sim->bnd.dat[seg->pid];
        struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

        fun_ctx.seg = si;
        fun_ctx.qud = seg->qud;

        double hxi = 0; // segment length
        int    nrm = asm_seg_dim(sim, seg, &hxi);

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                mx[i][j] = M[i][j] * hxi;

                for (int k = 0; k < 2; ++k) {
                    mnx[i][j][k] = MN[i][j][k] * hxi;
                }
            }
        }

        switch (cnd->type) {
            case CND_BND_DIR: {
                // dirichlet conditions
                if (log_add(&dir, si))
                    return -1;

                break;
            }
            case CND_BND_NEU: {
                // neumann conditions

                if (!ops.vneu)
                    continue;

                static double tta[2];
                struct val   *tta_v = &cnd->pps.neu.tta;

                // decompose parameters

                switch (tta_v->type) {
                    case VAL_NUM:
                        for (int k = 0; k < 2; ++k) {
                            tta[k] = tta_v->as.num;
                        }

                        break;
                    case VAL_FUN:
                        for (int k = 0; k < 2; ++k) {
                            fun_ctx.vtx = seg->vtx[k];
                            fun_ctx.ctx = tta_v->as.fun.ctx;
                            tta[k] = tta_v->as.fun.run(&fun_ctx, &(struct vec){
                                                                     .n = 2,
                                                                     .dat = vtx[fun_ctx.vtx].dat,
                                                                 });
                        }

                        break;
                    default:
                        return -ENOTSUP;
                }

                for (int i = 0; i < 2; ++i) {
                    int gi = seg->vtx[i]; // global i-index
                    int li = umsh_qud_vtx_loc(qud, gi);

                    int mui = MU[li];
                    int nui = NU[li];

                    double bi = 0;

                    for (int k = 0; k < 2; ++k) {
                        int lk = umsh_qud_vtx_loc(qud, seg->vtx[k]);
                        int muk = MU[lk];
                        int nuk = NU[lk];

                        double fix = 0;
                        double a = 0;
                        double b = 0;
                        double v = 0;

                        switch (nrm) {
                            case 1:                          // x-norm, x-fixed, vertical
                                v = vtx[seg->vtx[0]].dat[0]; // x
                                a = vtx[qud->vtx[0]].dat[0]; // x0
                                b = vtx[qud->vtx[1]].dat[0]; // x1

                                fix = F[muk](a, b, v) * F[mui](a, b, v);
                                bi += tta[k] * fix * mx[nuk][nui];

                                break;
                            case 2:                          // y-norm, y-fixed, horizontal
                                v = vtx[seg->vtx[0]].dat[1]; // y
                                a = vtx[qud->vtx[0]].dat[1]; // y0
                                b = vtx[qud->vtx[2]].dat[1]; // y1

                                fix = F[nuk](a, b, v) * F[nui](a, b, v);
                                bi += tta[k] * fix * mx[muk][mui];

                                break;
                        }
                    }

                    ops.vneu->dat[gi] += bi;
                }

                break;
            }
            case CND_BND_ROB: {
                // robin conditions

                if (ops.vrob == NULL && ops.mrob == NULL)
                    continue;

                static double bet[2];
                static double ext[2];

                struct val *bet_v = &cnd->pps.rob.bet;
                struct val *ext_v = &cnd->pps.rob.ext;

                // decompose parameters

                switch (bet_v->type) {
                    case VAL_NUM:
                        for (int k = 0; k < 2; ++k) {
                            bet[k] = bet_v->as.num;
                        }

                        break;
                    case VAL_FUN:
                        for (int k = 0; k < 2; ++k) {
                            fun_ctx.vtx = seg->vtx[k];
                            fun_ctx.ctx = bet_v->as.fun.ctx;
                            bet[k] = bet_v->as.fun.run(&fun_ctx, &(struct vec){
                                                                     .n = 2,
                                                                     .dat = vtx[fun_ctx.vtx].dat,
                                                                 });
                        }

                        break;
                    default:
                        return -ENOTSUP;
                }

                switch (ext_v->type) {
                    case VAL_NUM:
                        for (int k = 0; k < 2; ++k) {
                            ext[k] = ext_v->as.num;
                        }

                        break;
                    case VAL_FUN:
                        for (int k = 0; k < 2; ++k) {
                            fun_ctx.vtx = seg->vtx[k];
                            fun_ctx.ctx = ext_v->as.fun.ctx;
                            ext[k] = ext_v->as.fun.run(&fun_ctx, &(struct vec){
                                                                     .n = 2,
                                                                     .dat = vtx[fun_ctx.vtx].dat,
                                                                 });
                        }

                        break;
                    default:
                        return -ENOTSUP;
                }

                for (int i = 0; i < 2; ++i) {
                    int gi = seg->vtx[i]; // global i-index

                    int mui = MU[i];
                    // int nui = NU[i];

                    double bi = 0;

                    for (int j = 0; j < 2; ++j) {
                        int gj = seg->vtx[j]; // global j-index

                        int muj = MU[j];
                        // int nuj = NU[j];

                        double mij = 0;
                        double bj = 0;

                        for (int k = 0; k < 2; ++k) {
                            int muk = MU[k];
                            // int nuk = NU[k];

                            mij += bet[k] * mnx[muk][muj][mui];
                            bj += ext[k] * mnx[muk][muj][mui];
                        }

                        if (ops.mrob) {
                            mtx_inc(ops.mrob, gi, gj, mij);
                        }

                        bi += bet[j] * bj;
                    }

                    if (ops.vrob) {
                        ops.vrob->dat[gi] += bi;
                    }
                }

                break;
            }
        }
    }

    // apply dirichlet conditions
    if (ops.vdir || ops.mdir) {
        log_rst(&dir);

        for (int si = 0; !log_adv(&dir, &si);) {
            struct seg     *seg = &sim->msh->seg.dat[si];
            struct bnd     *bnd = &sim->bnd.dat[seg->pid];
            struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];
            struct val     *tgt_val = &cnd->pps.dir.tgt;

            fun_ctx.seg = si;
            fun_ctx.qud = seg->qud;

            double tgt_n = tgt_val->as.num;

            for (int i = 0; i < 2; ++i) {
                int gi = seg->vtx[i]; // global i-index

                if (ops.mdir) {
                    smtx_row_rst(ops.mdir, gi);
                    ops.mdir->dr[gi] = 1.0;
                }

                if (ops.vdir) {
                    fun_ctx.vtx = gi;

                    if (tgt_val->type == VAL_FUN) {
                        fun_ctx.ctx = tgt_val->as.fun.ctx;
                        ops.vdir->dat[gi] = tgt_val->as.fun.run(&fun_ctx, &(struct vec){
                                                                              .n = 2,
                                                                              .dat = vtx[gi].dat,
                                                                          });
                    } else {
                        ops.vdir->dat[gi] = tgt_n;
                    }
                }
            }
        }

        if (errno != ENOENT) {
            return -1;
        }

        errno = 0;
    }

    return 0;
}

static void twk_wgt(void* ctx, double hop, struct dif_ops* ops)
{
    struct sim_fun_ctx *fun_ctx = (struct sim_fun_ctx*)ctx;
    fun_ctx->sim->slv->run.wgt[0]->dat[ops->var] = hop;
}

int fem_std_lin_c2d_new(struct sim *sim, struct fem_std_ctx *ctx)
{
    struct v2d *vtx = sim->msh->vtx.v2d.dat;

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        struct qud *qud = &sim->msh->qud.dat[qi];
        struct obj *obj = &sim->obj.dat[qud->pid];
        struct mat *mat = &sim->mat.dat[obj->mat];

        if (mat->lam.type != VAL_FUN) {
            continue;
        }

        int v0 = qud->vtx[0];
        int v3 = qud->vtx[3];

        double hx = vtx[v3].dat[0] - vtx[v0].dat[0];
        double hy = vtx[v3].dat[1] - vtx[v0].dat[1];

        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    gnx[i][j][k] = GN[i][j][k] / hx;
                    gny[i][j][k] = GN[i][j][k] / hy;
                    mnx[i][j][k] = MN[i][j][k] * hx;
                    mny[i][j][k] = MN[i][j][k] * hy;
                }
            }

        double         dlam[4];
        struct dif_ops dif_ops = {0};

        for (int k = 0; k < 4; ++k) {
            int gk = qud->vtx[k];

            dif_ops.var = gk;
            dif_ops.hop = 0.000001;
            dif_ops.twk = twk_wgt;
            fun_ctx.vtx = gk;
            fun_ctx.qud = qi;
            fun_ctx.ctx = mat->lam.as.fun.ctx;

            struct vec vw = {.dat = vtx[gk].dat, .n = 2};
            dlam[k] = mat->lam.as.fun.dif(&fun_ctx, mat->lam.as.fun.run, &vw, &dif_ops);
        }

        for (int i = 0; i < 4; ++i) {
            int gi = qud->vtx[i];

            int mui = MU[i];
            int nui = NU[i];

            double bi = 0;

            for (int j = 0; j < 4; ++j) {
                int gj = qud->vtx[j];

                int muj = MU[j];
                int nuj = NU[j];

                double mij = 0;
                double bij = 0;

                for (int k = 0; k < 4; ++k) {
                    int gk = qud->vtx[k];

                    int muk = MU[k];
                    int nuk = NU[k];

                    double gln = gnx[muk][muj][mui] *
                                 mny[nuk][nuj][nui] +
                                 mnx[muk][muj][mui] *
                                 gny[nuk][nuj][nui];

                    mij += ctx->w0.dat[gk] * dlam[j] * gln;
                    bij += ctx->w0.dat[gk] * dlam[k] * gln;
                }

                bi += bij * ctx->w0.dat[gj];

                mtx_inc(&ctx->mtx, gi, gj, mij * DPH);
            }

            ctx->vec.dat[gi] += bi * DPH;
        }
    }

    return 0;
}
