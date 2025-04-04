#include <math.h>

#include <numx/com/log.h>

#include "fem.h"

#define C 10e10

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

static const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
static const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
static const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

static int asm_hxd(struct sim *sim, struct vec *wgt, int h);

static int asm_qud_dir(struct sim *sim, struct vec *wgt, int q);
static int asm_qud_neu(struct sim *sim, struct vec *wgt, int q);
static int asm_qud_rob(struct sim *sim, struct vec *wgt, int q);

static int asm_mov_mtx(struct smtx *m, int i, int j, double v);

int fem_ell_lin_std_slv(struct sim *sim, struct vec *wgt)
{
    int hl = sim->msh->hxd.len;
    int ql = sim->msh->qud.len;

    for (int i = 0; i < hl; ++i)
        if (asm_hxd(sim, wgt, i))
            return -1;

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    struct qud     *q = sim->msh->qud.dat;
    struct bnd     *b = sim->bnd.dat;
    struct cnd_bnd *c = sim->cnd_bnd.dat;

    for (int i = 0; i < ql; ++i) {
        struct qud     *qud = &q[i];
        struct bnd     *bnd = &b[qud->pid];
        struct cnd_bnd *cnd = &c[bnd->cnd];

        switch (cnd->type) {
            case CND_BND_DIR:
                if (log_add(&dir, i)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_NEU:
                if (asm_qud_neu(sim, wgt, i)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_ROB:
                if (asm_qud_rob(sim, wgt, i)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
        }
    }

    log_rst(&dir);

    for (int i = 0; !log_adv(&dir, &i);)
        if (asm_qud_dir(sim, wgt, i)) {
            log_cls(&dir);
            return -1;
        }

    log_cls(&dir);

    if (errno != ENOENT)
        return -1;

    errno = 0;

    switch (sim->fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&sim->fem->mtx, wgt, &sim->fem->vec, sim->fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int asm_hxd(struct sim *sim, struct vec *wgt, int h)
{
    static double lam[8];
    static double gam[8];
    static double src[8];

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

    struct hxd *hxd = &sim->msh->hxd.dat[h];
    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];
    struct val *val = &sim->src.dat[obj->src];

    struct vtx *vtx = sim->msh->vtx.dat;

    if (mat->lam.type == VAL_FUN) {
        fun lf = mat->lam.as.fun;

        for (int k = 0; k < 8; ++k)
            lam[k] = lf(sim, wgt, hxd->vtx[k]);
    }

    if (mat->gam.type == VAL_FUN) {
        fun gf = mat->gam.as.fun;

        for (int k = 0; k < 8; ++k)
            gam[k] = gf(sim, wgt, hxd->vtx[k]);
    }

    if (val->type == VAL_FUN) {
        fun sf = val->as.fun;

        for (int k = 0; k < 8; ++k)
            src[k] = sf(sim, wgt, hxd->vtx[k]);
    }

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
                mij += mat->lam.as.num * (gx[muj][mui] *
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
                mij += mat->gam.as.num * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);

            asm_mov_mtx(&sim->fem->mtx, gi, gj, mij);
        }

        double bi = 0;

        if (val->type == VAL_FUN)
            for (int k = 0; k < 8; ++k) {
                int muk = MU[k];
                int nuk = NU[k];
                int ttk = TT[k];

                bi += src[k] * (mx[muk][mui] * my[nuk][nui] * mz[ttk][tti]);
            }
        else
            bi = val->as.num * hx * hy * hz / 8;

        sim->fem->vec.dat[gi] += bi;
    }

    return 0;
}

static int asm_qud_dir(struct sim *sim, struct vec *wgt, int q)
{
    struct qud     *qud = &sim->msh->qud.dat[q];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[sim->bnd.dat[qud->pid].cnd];

    struct smtx *mtx = &sim->fem->mtx;
    struct vec  *vec = &sim->fem->vec;

    fun    ftgt = cnd->pps.dir.tgt.as.fun;
    double ntgt = cnd->pps.dir.tgt.as.num;

    for (int i = 0; i < 4; ++i) {
        int gi = qud->vtx[i];

        mtx->dr[gi] = C;

        if (cnd->pps.dir.tgt.type == VAL_FUN)
            vec->dat[gi] = C * ftgt(sim, wgt, gi);
        else
            vec->dat[gi] = C * ntgt;
    }

    return 0;
}

static int asm_qud_neu(struct sim *sim, struct vec *wgt, int q)
{
    static double tta[4];
    static double mxi[2][2];
    static double mzt[2][2];

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[sim->bnd.dat[qud->pid].cnd];

    fun    ftta = cnd->pps.neu.tta.as.fun;
    double ntta = cnd->pps.neu.tta.as.num;

    if (cnd->pps.neu.tta.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            tta[i] = ftta(sim, wgt, qud->vtx[i]);

    int    nrm = msh_qud_nrm(sim->msh, q);
    double hxi = 0;
    double hzt = 0;

    struct vtx *vtx = sim->msh->vtx.dat;

    switch (nrm) {
        case 0:
            hxi = fabs(vtx[qud->vtx[0]].y - vtx[qud->vtx[2]].y);
            hzt = fabs(vtx[qud->vtx[0]].z - vtx[qud->vtx[2]].z);

            break;
        case 1:
            hxi = fabs(vtx[qud->vtx[0]].x - vtx[qud->vtx[2]].x);
            hzt = fabs(vtx[qud->vtx[0]].z - vtx[qud->vtx[2]].z);

            break;
        case 2:
            hxi = fabs(vtx[qud->vtx[0]].x - vtx[qud->vtx[2]].x);
            hzt = fabs(vtx[qud->vtx[0]].y - vtx[qud->vtx[2]].y);

            break;
    }

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

        if (cnd->pps.neu.tta.type == VAL_FUN)
            for (int k = 0; k < 4; ++k) {
                int muk = MU[k];
                int nuk = NU[k];

                bi += tta[k] * mxi[muk][mui] * mzt[nuk][nui];
            }
        else
            bi = ntta * hxi * hzt / 4;

        sim->fem->vec.dat[gi] += bi;
    }

    return 0;
}

static int asm_qud_rob(struct sim *sim, struct vec *wgt, int q)
{
    static double bet[4];
    static double ext[4];

    static double mxi[2][2];
    static double mzt[2][2];

    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    struct qud     *qud = &sim->msh->qud.dat[q];
    struct cnd_bnd *cnd = &sim->cnd_bnd.dat[sim->bnd.dat[qud->pid].cnd];

    fun    fbet = cnd->pps.rob.bet.as.fun;
    double nbet = cnd->pps.rob.bet.as.num;

    fun    fext = cnd->pps.rob.ext.as.fun;
    double next = cnd->pps.rob.ext.as.num;

    if (cnd->pps.rob.bet.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            bet[i] = fbet(sim, wgt, qud->vtx[i]);

    if (cnd->pps.rob.ext.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            ext[i] = fext(sim, wgt, qud->vtx[i]);

    int    nrm = msh_qud_nrm(sim->msh, q);
    double hxi = 0;
    double hzt = 0;

    struct vtx *vtx = sim->msh->vtx.dat;

    int v0 = qud->vtx[0];
    int v2 = qud->vtx[2];

    switch (nrm) {
        case 0:
            hxi = fabs(vtx[v0].y - vtx[v2].y);
            hzt = fabs(vtx[v0].z - vtx[v2].z);

            break;
        case 1:
            hxi = fabs(vtx[v0].x - vtx[v2].x);
            hzt = fabs(vtx[v0].z - vtx[v2].z);

            break;
        case 2:
            hxi = fabs(vtx[v0].x - vtx[v2].x);
            hzt = fabs(vtx[v0].y - vtx[v2].y);

            break;
    }

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            mxi[i][j] = M[i][j] * hxi;
            mzt[i][j] = M[i][j] * hzt;

            for (int k = 0; k < 2; ++k) {
                mnxi[i][j][k] = MN[i][j][k] * hxi;
                mnzt[i][j][k] = MN[i][j][k] * hzt;
            }
        }

    struct smtx *mtx = &sim->fem->mtx;
    struct vec  *vec = &sim->fem->vec;

    switch (cnd->pps.rob.bet.type) {
        case VAL_NUM:
            for (int i = 0; i < 4; ++i) {
                int gi = qud->vtx[i];
                int mui = MU[i];
                int nui = NU[i];

                for (int j = 0; j < 4; ++j) {
                    int gj = qud->vtx[j];
                    int muj = MU[j];
                    int nuj = NU[j];

                    double mij = nbet * (mxi[muj][mui] * mzt[nuj][nui]);

                    if (asm_mov_mtx(mtx, gi, gj, mij))
                        return -1;
                }
            }

            switch (cnd->pps.rob.ext.type) {
                case VAL_NUM:
                    for (int i = 0; i < 4; ++i)
                        vec->dat[qud->vtx[i]] += nbet * next * hxi * hzt / 4;

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

                        vec->dat[qud->vtx[i]] += nbet * bi;
                    }

                    break;
            }

            break;
        case VAL_FUN:
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

                        mij += bet[k] * (mnxi[muk][muj][mui] * mnzt[nuk][nuj][nui]);
                    }

                    if (asm_mov_mtx(mtx, gi, gj, mij))
                        return -1;
                }
            }

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

                        vec->dat[qud->vtx[i]] += next * bi;
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

static int asm_mov_mtx(struct smtx *m, int i, int j, double v)
{
    if (i < j) {
        int p = m->ia[j];

        while (m->ja[p] < i)
            ++p;

        m->ur[p] += v;

        return 0;
    }

    if (j < i) {
        int p = m->ia[i];

        while (m->ja[p] < j)
            ++p;

        m->lr[p] += v;

        return 0;
    }

    m->dr[i] += v;

    return 0;
}
