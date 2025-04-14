#include <assert.h>
#include <math.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/sim.h>

#include "fem_lin.h"

static int fem_lin_mtx_lam_hxd_asm(struct sim *sim, struct smtx *mtx, int h);
static int fem_lin_mtx_gam_hxd_asm(struct sim *sim, struct smtx *mtx, int h);
static int fem_lin_mtx_sig_hxd_asm(struct sim *sim, struct smtx *mtx, int h);
static int fem_lin_mtx_chi_hxd_asm(struct sim *sim, struct smtx *mtx, int h);

static int fem_lin_mtx_rob_qud_asm(struct sim *sim, struct smtx *mtx, int q);
static int fem_lin_mtx_dir_qud_asm(struct sim *sim, struct smtx *mtx, int q);

static int fem_lin_vec_src_hxd_asm(struct sim *sim, struct vec *vec, int h);
static int fem_lin_vec_dir_qud_asm(struct sim *sim, struct vec *vec, int q);
static int fem_lin_vec_neu_qud_asm(struct sim *sim, struct vec *vec, int q);
static int fem_lin_vec_rob_qud_asm(struct sim *sim, struct vec *vec, int q);

int fem_lin_ell_asm(struct sim *sim)
{
    assert(sim);

    struct fem *fem = (struct fem *)sim->slv;

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
        if (fem_lin_mtx_lam_hxd_asm(sim, &fem->prv.ell, hi))
            return -1;

        if (fem_lin_mtx_gam_hxd_asm(sim, &fem->prv.ell, hi))
            return -1;

        if (fem_lin_vec_src_hxd_asm(sim, &fem->prv.vec, hi))
            return -1;
    }

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        switch (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type) {
            case CND_BND_DIR:
                if (log_add(&dir, qi))
                    return -1;

                break;
            case CND_BND_NEU:
                if (fem_lin_vec_neu_qud_asm(sim, &fem->prv.vec, qi))
                    return -1;

                break;
            case CND_BND_ROB:
                if (fem_lin_mtx_rob_qud_asm(sim, &fem->prv.ell, qi))
                    return -1;

                if (fem_lin_vec_rob_qud_asm(sim, &fem->prv.vec, qi))
                    return -1;

                break;
        }

    log_rst(&dir);

    for (int qi = 0; !log_adv(&dir, &qi);) {
        if (fem_lin_mtx_dir_qud_asm(sim, &fem->prv.ell, qi))
            return -1;

        if (fem_lin_vec_dir_qud_asm(sim, &fem->prv.vec, qi))
            return -1;
    }

    if (errno != ENOENT)
        return -1;

    errno = 0;
    return 0;
}

int fem_lin_mtx_lam_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (fem_lin_mtx_lam_hxd_asm(sim, mtx, hi))
            return -1;

    return 0;
}

int fem_lin_mtx_gam_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (fem_lin_mtx_gam_hxd_asm(sim, mtx, hi))
            return -1;

    return 0;
}

int fem_lin_mtx_sig_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (fem_lin_mtx_sig_hxd_asm(sim, mtx, hi))
            return -1;

    return 0;
}

int fem_lin_mtx_chi_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (fem_lin_mtx_chi_hxd_asm(sim, mtx, hi))
            return -1;

    return 0;
}

int fem_lin_mtx_rob_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        if (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type == CND_BND_ROB)
            if (fem_lin_mtx_rob_qud_asm(sim, mtx, qi))
                return -1;

    return 0;
}

int fem_lin_mtx_dir_asm(struct sim *sim, struct smtx *mtx)
{
    assert(sim);
    assert(mtx);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        if (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type == CND_BND_DIR)
            if (fem_lin_mtx_dir_qud_asm(sim, mtx, qi))
                return -1;

    return 0;
}

int fem_lin_vec_src_asm(struct sim *sim, struct vec *vec)
{
    assert(sim);
    assert(vec);

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi)
        if (fem_lin_vec_src_hxd_asm(sim, vec, hi))
            return -1;

    return 0;
}

int fem_lin_vec_neu_asm(struct sim *sim, struct vec *vec)
{
    assert(sim);
    assert(vec);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        if (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type == CND_BND_NEU)
            if (fem_lin_vec_neu_qud_asm(sim, vec, qi))
                return -1;

    return 0;
}

int fem_lin_vec_rob_asm(struct sim *sim, struct vec *vec)
{
    assert(sim);
    assert(vec);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        if (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type == CND_BND_ROB)
            if (fem_lin_vec_rob_qud_asm(sim, vec, qi))
                return -1;

    return 0;
}

int fem_lin_vec_dir_asm(struct sim *sim, struct vec *vec)
{
    assert(sim);
    assert(vec);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi)
        if (sim->cnd_bnd.dat[sim->bnd.dat[sim->msh->qud.dat[qi].pid].cnd].type == CND_BND_DIR)
            if (fem_lin_vec_dir_qud_asm(sim, vec, qi))
                return -1;

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

static int fem_lin_mtx_lam_hxd_asm(struct sim *sim, struct smtx *mtx, int h)
{
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

    static double lam[8];

    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    double lam_n = mat->lam.as.num;

    if (mat->lam.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            lam[k] = mat->lam.as.fun(sim, hxd->vtx[k]);

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

            if (mat->lam.type == VAL_FUN) {
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
            } else {
                mij += lam_n * (gx[muj][mui] *
                                   my[nuj][nui] *
                                   mz[ttj][tti] +
                                   mx[muj][mui] *
                                   gy[nuj][nui] *
                                   mz[ttj][tti] +
                                   mx[muj][mui] *
                                   my[nuj][nui] *
                                   gz[ttj][tti]);
            }

            if (smtx_add(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int fem_lin_mtx_gam_hxd_asm(struct sim *sim, struct smtx *mtx, int h)
{
    static double mx[2][2];
    static double my[2][2];
    static double mz[2][2];

    static double mnx[2][2][2];
    static double mny[2][2][2];
    static double mnz[2][2][2];

    static double gam[8];

    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    double gam_n = mat->gam.as.num;

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

                    mij += gam[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                }
            } else {
                mij += gam_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
            }

            if (smtx_add(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int fem_lin_mtx_sig_hxd_asm(struct sim *sim, struct smtx *mtx, int h)
{
    static double mx[2][2];
    static double my[2][2];
    static double mz[2][2];

    static double mnx[2][2][2];
    static double mny[2][2][2];
    static double mnz[2][2][2];

    static double sig[8];

    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    double sig_n = mat->sig.as.num;

    if (mat->sig.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            sig[k] = mat->sig.as.fun(sim, hxd->vtx[k]);

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

                    mij += sig[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                }
            } else {
                mij += sig_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
            }

            if (smtx_add(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int fem_lin_mtx_chi_hxd_asm(struct sim *sim, struct smtx *mtx, int h)
{
    static double mx[2][2];
    static double my[2][2];
    static double mz[2][2];

    static double mnx[2][2][2];
    static double mny[2][2][2];
    static double mnz[2][2][2];

    static double chi[8];

    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];

    double chi_n = mat->chi.as.num;

    if (mat->chi.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            chi[k] = mat->chi.as.fun(sim, hxd->vtx[k]);

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

                    mij += chi[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                }
            } else {
                mij += chi_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);
            }

            if (smtx_add(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int fem_lin_mtx_rob_qud_asm(struct sim *sim, struct smtx *mtx, int q)
{
    assert(sim);
    assert(mtx);

    static double mxi[2][2];
    static double mzt[2][2];
    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    static double bet[4];

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct bnd     *bnd = &sim->bnd.dat[qud->pid];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

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

            if (smtx_add(mtx, gi, gj, mij))
                return -1;
        }
    }

    return 0;
}

static int fem_lin_mtx_dir_qud_asm(struct sim *sim, struct smtx *mtx, int q)
{
    assert(sim);
    assert(mtx);

    struct qud *qud = &sim->msh->qud.dat[q];

    for (int i = 0; i < 4; ++i)
        mtx->dr[qud->vtx[i]] = C;

    return 0;
}

static int fem_lin_vec_src_hxd_asm(struct sim *sim, struct vec *vec, int h)
{
    assert(sim);
    assert(vec);

    static double mx[2][2];
    static double my[2][2];
    static double mz[2][2];

    static double src[8];

    struct vtx *vtx = sim->msh->vtx.dat;
    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct val *val = &sim->src.dat[obj->src];

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

static int fem_lin_vec_dir_qud_asm(struct sim *sim, struct vec *vec, int q)
{
    assert(sim);
    assert(vec);

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct bnd     *bnd = &sim->bnd.dat[qud->pid];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

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

static int fem_lin_vec_neu_qud_asm(struct sim *sim, struct vec *vec, int q)
{
    static double mxi[2][2];
    static double mzt[2][2];

    static double tta[4];

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct bnd     *bnd = &sim->bnd.dat[qud->pid];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

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
        } else {
            bi = tta_n * hxi * hzt / 4;
        }

        vec->dat[gi] += bi;
    }

    return 0;
}

static int fem_lin_vec_rob_qud_asm(struct sim *sim, struct vec *vec, int q)
{
    assert(sim);
    assert(vec);

    static double mxi[2][2];
    static double mzt[2][2];
    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    static double bet[4];
    static double ext[4];

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct bnd     *bnd = &sim->bnd.dat[qud->pid];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

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
