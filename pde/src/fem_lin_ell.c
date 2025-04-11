#include <assert.h>
#include <math.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/sim.h>

#include "fem_lin.h"

static int asm_hxd_vec(struct sim *, int);
static int asm_hxd_mtx(struct sim *, int);

static int asm_qud_dir_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_dir_mtx(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_neu_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_rob_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_rob_mtx(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd);

int fem_lin_ell_mtx_asm(struct sim *sim, struct vec *wgt)
{
    assert(sim);
    assert(wgt);

    sim->rt.wgt = wgt;

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (asm_hxd_mtx(sim, hi))
            return -1;

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    struct qud     *qud = sim->msh->qud.dat;
    struct bnd     *bnd = sim->bnd.dat;
    struct cnd_bnd *cnd = sim->cnd_bnd.dat;

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        struct qud     *q = &qud[qi];
        struct bnd     *b = &bnd[q->pid];
        struct cnd_bnd *c = &cnd[b->cnd];

        switch (c->type) {
            case CND_BND_DIR:
                if (log_add(&dir, qi)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_NEU:
                break;
            case CND_BND_ROB:
                if (asm_qud_rob_mtx(sim, q, c)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
        }
    }

    log_rst(&dir);

    for (int qi = 0; !log_adv(&dir, &qi);) {
        struct qud     *q = &qud[qi];
        struct bnd     *b = &bnd[q->pid];
        struct cnd_bnd *c = &cnd[b->cnd];

        if (asm_qud_dir_mtx(sim, q, c)) {
            log_cls(&dir);
            return -1;
        }
    }

    log_cls(&dir);

    if (errno != ENOENT)
        return -1;

    errno = 0;
    return 0;
}

int fem_lin_ell_vec_asm(struct sim *sim, struct vec *wgt)
{
    assert(sim);
    assert(wgt);

    sim->rt.wgt = wgt;

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (asm_hxd_vec(sim, hi))
            return -1;

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    struct qud     *qud = sim->msh->qud.dat;
    struct bnd     *bnd = sim->bnd.dat;
    struct cnd_bnd *cnd = sim->cnd_bnd.dat;

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        struct qud     *q = &qud[qi];
        struct bnd     *b = &bnd[q->pid];
        struct cnd_bnd *c = &cnd[b->cnd];

        switch (c->type) {
            case CND_BND_DIR:
                if (log_add(&dir, qi)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_NEU:
                if (asm_qud_neu_vec(sim, q, c)) {
                    log_cls(&dir);
                    return -1;
                }
                break;
            case CND_BND_ROB:
                if (asm_qud_rob_vec(sim, q, c)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
        }
    }

    log_rst(&dir);

    for (int qi = 0; !log_adv(&dir, &qi);) {
        struct qud     *q = &qud[qi];
        struct bnd     *b = &bnd[q->pid];
        struct cnd_bnd *c = &cnd[b->cnd];

        if (asm_qud_dir_vec(sim, q, c)) {
            log_cls(&dir);
            return -1;
        }
    }

    log_cls(&dir);

    if (errno != ENOENT)
        return -1;

    errno = 0;
    return 0;
}

static int asm_hxd_vec(struct sim *sim, int hi)
{
    assert(sim);

    struct fem *fem = (struct fem *)sim->slv;
    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[hi];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct val *val = &sim->src.dat[obj->src];

    static double src[8];

    static double mx[2][2];
    static double my[2][2];
    static double mz[2][2];

    double src_n = val->as.num;

    if (val->type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            src[k] = val->as.fun(sim, hxd->vtx[k]);

    int v0 = hxd->vtx[0];
    int v7 = hxd->vtx[7];

    double hx = vtx[v7].x - vtx[v0].x;
    double hy = vtx[v7].y - vtx[v0].y;
    double hz = vtx[v7].z - vtx[v0].z;

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

        if (val->type == VAL_FUN)
            for (int k = 0; k < 8; ++k) {
                int muk = MU[k];
                int nuk = NU[k];
                int ttk = TT[k];

                bi += src[k] * (mx[muk][mui] * my[nuk][nui] * mz[ttk][tti]);
            }
        else
            bi = src_n * hx * hy * hz / 8;

        fem->prv.vec.dat[gi] += bi;
    }

    return 0;
}

static int asm_hxd_mtx(struct sim *sim, int hi)
{
    assert(sim);

    struct fem *fem = (struct fem *)sim->slv;
    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[hi];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    static double lam[8];
    static double gam[8];

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

    double lam_n = mat->lam.as.num;
    double gam_n = mat->gam.as.num;

    if (mat->lam.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            lam[k] = mat->lam.as.fun(sim, hxd->vtx[k]);

    if (mat->gam.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            gam[k] = mat->gam.as.fun(sim, hxd->vtx[k]);

    int v0 = hxd->vtx[0];
    int v7 = hxd->vtx[7];

    double hx = vtx[v7].x - vtx[v0].x;
    double hy = vtx[v7].y - vtx[v0].y;
    double hz = vtx[v7].z - vtx[v0].z;

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

            if (mat->lam.type == VAL_FUN)
                for (int k = 0; k < 8; ++k) {
                    int muk = MU[k];
                    int nuk = NU[k];
                    int ttk = TT[k];

                    mij += lam[k] * (gnx[muj][mui] *
                                        mny[nuk][nuj][nui] *
                                        mnz[ttk][ttj][tti] +
                                        mnx[muk][muj][mui] *
                                        gny[nuj][nui] *
                                        mnz[ttk][ttj][tti] +
                                        mnx[muk][muj][mui] *
                                        mny[nuk][nuj][nui] *
                                        gnz[ttj][tti]);
                }
            else
                mij += lam_n * (gx[muj][mui] *
                                   my[nuj][nui] *
                                   mz[ttj][tti] +
                                   mx[muj][mui] *
                                   gy[nuj][nui] *
                                   mz[ttj][tti] +
                                   mx[muj][mui] *
                                   my[nuj][nui] *
                                   gz[ttj][tti]);

            if (mat->gam.type == VAL_FUN)
                for (int k = 0; k < 8; ++k) {
                    int muk = MU[k];
                    int nuk = NU[k];
                    int ttk = TT[k];

                    mij += gam[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                }
            else
                mij += gam_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);

            asm_mov_mtx(&fem->prv.ell, gi, gj, mij);
        }
    }

    return 0;
}

static int asm_qud_dir_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd)
{
    assert(sim);

    struct fem *fem = (struct fem *)sim->slv;
    struct vec *vec = &fem->prv.vec;

    double tgt_n = cnd->pps.dir.tgt.as.num;

    for (int i = 0; i < 4; ++i) {
        int gi = qud->vtx[i];

        if (cnd->pps.dir.tgt.type == VAL_FUN)
            vec->dat[gi] = C * cnd->pps.dir.tgt.as.fun(sim, gi);
        else
            vec->dat[gi] = C * tgt_n;
    }

    return 0;
}

static int asm_qud_dir_mtx(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd)
{
    assert(sim);
    assert(cnd);

    struct fem  *fem = (struct fem *)sim->slv;
    struct smtx *mtx = &fem->prv.ell;

    for (int i = 0; i < 4; ++i)
        mtx->dr[qud->vtx[i]] = C;

    return 0;
}

static int asm_qud_dim(struct sim *sim, struct qud *qud, double *hxi, double *hzt);

static int asm_qud_neu_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd)
{
    static double tta[4];
    static double mxi[2][2];
    static double mzt[2][2];

    struct fem *fem = (struct fem *)sim->slv;
    struct vec *vec = &fem->prv.vec;

    double tta_n = cnd->pps.neu.tta.as.num;

    if (cnd->pps.neu.tta.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            tta[i] = cnd->pps.neu.tta.as.fun(sim, qud->vtx[i]);

    double hxi = 0;
    double hzt = 0;

    if (asm_qud_dim(sim, qud, &hxi, &hzt))
        return -1;

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            mxi[i][j] = M[i][j] * hxi;
            mzt[i][j] = M[i][j] * hzt;
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

                bi += tta[k] * mxi[muk][mui] * mzt[nuk][nui];
            }
        } else
            bi = tta_n * hxi * hzt / 4;

        vec->dat[gi] += bi;
    }

    return 0;
}

static int asm_qud_rob_vec(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd)
{
    assert(sim);
    assert(qud);
    assert(cnd);

    static double bet[4];
    static double ext[4];

    static double mxi[2][2];
    static double mzt[2][2];

    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    struct fem *fem = (struct fem *)sim->slv;
    struct vec *vec = &fem->prv.vec;

    double bet_n = cnd->pps.rob.bet.as.num;
    double ext_n = cnd->pps.rob.ext.as.num;

    if (cnd->pps.rob.bet.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            bet[i] = cnd->pps.rob.bet.as.fun(sim, qud->vtx[i]);

    if (cnd->pps.rob.ext.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            ext[i] = cnd->pps.rob.ext.as.fun(sim, qud->vtx[i]);

    double hxi = 0;
    double hzt = 0;

    if (asm_qud_dim(sim, qud, &hxi, &hzt))
        return -1;

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            mxi[i][j] = M[i][j] * hxi;
            mzt[i][j] = M[i][j] * hzt;

            for (int k = 0; k < 2; ++k) {
                mnxi[i][j][k] = MN[i][j][k] * hxi;
                mnzt[i][j][k] = MN[i][j][k] * hzt;
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

                            bi += ext[k] * mxi[muk][mui] * mzt[nuk][nui];
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

                            bi += bet[k] * mxi[muk][mui] * mzt[nuk][nui];
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

                                bi += bet[k] * ext[j] * mnxi[muk][muj][mui] * mnzt[nuk][nuj][nui];
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

static int asm_qud_rob_mtx(struct sim *sim, struct qud *qud, struct cnd_bnd *cnd)
{
    assert(sim);
    assert(qud);
    assert(cnd);

    static double bet[4];

    static double mxi[2][2];
    static double mzt[2][2];
    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    struct fem  *fem = (struct fem *)sim->slv;
    struct smtx *mtx = &fem->prv.ell;

    double bet_n = cnd->pps.rob.bet.as.num;

    if (cnd->pps.rob.bet.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            bet[i] = cnd->pps.rob.bet.as.fun(sim, qud->vtx[i]);

    double hxi = 0;
    double hzt = 0;

    if (asm_qud_dim(sim, qud, &hxi, &hzt))
        return -1;

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            mxi[i][j] = M[i][j] * hxi;
            mzt[i][j] = M[i][j] * hzt;

            for (int k = 0; k < 2; ++k) {
                mnxi[i][j][k] = MN[i][j][k] * hxi;
                mnzt[i][j][k] = MN[i][j][k] * hzt;
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
                    mij = bet_n * (mxi[muj][mui] * mzt[nuj][nui]);
                    break;
                case VAL_FUN:
                    for (int k = 0; k < 4; ++k) {
                        int muk = MU[k];
                        int nuk = NU[k];

                        mij += bet[k] * (mnxi[muk][muj][mui] * mnzt[nuk][nuj][nui]);
                    }

                    break;
            }

            if (asm_mov_mtx(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int asm_qud_dim(struct sim *sim, struct qud *qud, double *hxi, double *hzt)
{
    double dat[3];

    struct vec nrm = {
        .n = 3,
        .dat = dat,
    };

    if (msh_qud_nrm(sim->msh, qud, &nrm))
        return -1;

    struct vtx *vtx = sim->msh->vtx.dat;

    if (nrm.dat[0] != 0) {
        *hxi = fabs(vtx[qud->vtx[0]].y - vtx[qud->vtx[3]].y);
        *hzt = fabs(vtx[qud->vtx[0]].z - vtx[qud->vtx[3]].z);
    } else if (nrm.dat[1] != 0) {
        *hxi = fabs(vtx[qud->vtx[0]].x - vtx[qud->vtx[3]].x);
        *hzt = fabs(vtx[qud->vtx[0]].z - vtx[qud->vtx[3]].z);
    } else {
        *hxi = fabs(vtx[qud->vtx[0]].x - vtx[qud->vtx[3]].x);
        *hzt = fabs(vtx[qud->vtx[0]].y - vtx[qud->vtx[3]].y);
    }

    return 0;
}

static int slv_std(struct sim *, struct vec *);
static int slv_non(struct sim *, struct vec *);

int fem_lin_ell_slv(struct sim *sim)
{
    assert(sim);

    int r = 0;

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  wgt;

    if ((r = vec_new(&wgt, fem->prv.vec.n)))
        goto end;

    sim->rt.tm = 0;

    switch (fem->ops.bss) {
        case FEM_BSS_LIN:
            sim->slv->apx.ctx = 0;
            sim->slv->apx.run = fem_lin_apx;

            switch (fem->ops.mod) {
                case FEM_STD:
                    if ((r = slv_std(sim, &wgt)))
                        goto end;

                    break;
                case FEM_NON:
                    if ((r = slv_non(sim, &wgt)))
                        goto end;

                    break;
                case FEM_HMC:
                    errno = ENOTSUP;
                    return -1;
            }

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    if (fem->slv.itr.run)
        fem->slv.itr.run(fem->slv.itr.ctx, sim);

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim);

end:
    vec_cls(&wgt);

    return r;
}

static int slv_std(struct sim *sim, struct vec *wgt)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem_lin_ell_mtx_asm(sim, wgt))
        return -1;

    if (fem_lin_ell_vec_asm(sim, wgt))
        return -1;

    switch (fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.ell, wgt, &fem->prv.vec, &fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int slv_non(struct sim *sim, struct vec *wgt)
{
    int r = 0;

    struct fem  *fem = (struct fem *)sim->slv;
    struct smtx *mtx = &fem->prv.ell;
    struct vec  *vec = &fem->prv.vec;
    struct vec   tmp;

    if ((r = vec_new(&tmp, vec->n)))
        goto end;

    int    max = fem->ops.non.max;
    double err = fem->ops.non.err;
    double cur = 1;
    double tn = 0;
    double bn = 0;

    for (int i = 0; i < max && cur > err; ++i) {
        if (i > 0)
            switch (fem->ops.iss.mod) {
                case ISS_BCG:
                    if ((r = iss_bcg_slv(mtx, wgt, vec, &fem->ops.iss.ops.bcg)))
                        goto end;

                    break;
                default:
                    errno = ENOTSUP;
                    r = -1;
                    goto end;
            }

        if ((r = fem_lin_ell_mtx_asm(sim, wgt)))
            goto end;

        if ((r = fem_lin_ell_vec_asm(sim, wgt)))
            goto end;

        mtx_vmlt(mtx, wgt, &tmp);
        vec_cmb(&tmp, vec, &tmp, -1);
        vec_nrm(&tmp, &tn);
        vec_nrm(vec, &bn);

        cur = tn / bn;

        fem->ops.non.rt.itr = i;
        fem->ops.non.rt.err = cur;

        if (fem->ops.non.itr.run)
            fem->ops.non.itr.run(fem->ops.non.itr.ctx, sim);
    }

end:
    vec_cls(&tmp);

    return r;
}
