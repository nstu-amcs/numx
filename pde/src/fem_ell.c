#include <assert.h>
#include <stdlib.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>

#define C 10e10

const double G[2][2] = {
    {1.0,  -1.0},
    {-1.0, 1.0 }
};

const double M[2][2] = {
    {2.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 6.0}
};

const double GN[2][2] = {
    {1.0 / 2.0,  -1.0 / 2.0},
    {-1.0 / 2.0, 1.0 / 2.0 }
};

const double MN[2][2][2] = {
    {{1.0 / 4.0, 1.0 / 12.0},  {1.0 / 12.0, 1.0 / 12.0}},
    {{1.0 / 12.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 4.0} }
};

const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

const int DM[8] = {0, 1, 3, 2, 4, 5, 7, 6};

static double fabs(double a)
{
    return a < 0 ? -a : a;
}

static int pfl(struct sim *sim, struct smtx *m);
static int slv(struct sim *sim, struct smtx *m, struct vec *b, struct vec *q);

int fem_ell_slv(struct sim *sim, struct vec *q)
{
    assert(sim);
    assert(q);

    int r = 0;

    int    max = sim->ops.ell.ops.non.max;
    double err = sim->ops.ell.ops.non.err;
    double cur = 1;

    struct smtx m;
    struct vec  b;
    struct vec  t;

    double tn = 0;
    double bn = 0;

    if ((r = pfl(sim, &m)))
        goto end;

    if ((r = vec_new(q, m.pps.n)))
        goto end;

    if ((r = vec_new(&t, m.pps.n)))
        goto end;

    if ((r = vec_new(&b, m.pps.n)))
        goto end;

    for (int i = 0; i < max && cur > err; ++i) {
        if ((r = slv(sim, &m, &b, q)))
            goto end;

        mtx_vmlt(&m, q, &t);

        vec_cmb(&t, &b, &t, -1);
        vec_nrm(&t, &tn);
        vec_nrm(&b, &bn);

        cur = tn / bn;
    }

end:
    mtx_cls(&m);
    vec_cls(&b);
    vec_cls(&t);

    return r;
}

static int pfl(struct sim *sim, struct smtx *m)
{
    int n = sim->msh->vtx.len;
    int z = 0;

    struct ilog *map = malloc(sizeof(struct ilog) * n);

    if (!map)
        return -1;

    for (int i = 0; i < n; ++i) {
        log_new(&map[i]);

        map[i].dup = false;
        map[i].srt = true;
        map[i].cmp.run = iasc;
    }

    for (int i = 0; i < sim->msh->hxd.len; ++i) {
        int *vtx = sim->msh->hxd.dat[i].vtx;

        for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 8; ++k)
                if (vtx[k] < vtx[j] && !log_add(&map[vtx[j]], vtx[k]))
                    z += 1;
    }

    if (mtx_new(m, ((struct smtx_pps){n, z}))) {
        free(map);
        return -1;
    }

    for (int i = 0, e = 0; i < n; ++i) {
        m->ia[i] = e;

        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e++)
            m->ja[e] = j;

        log_cls(&map[i]);
    }

    m->ia[n] = z;

    free(map);

    return 0;
}

static int asm_hxd(struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct hxd *hxd);

static int asm_qud_dir(
    struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_neu(struct sim *sim, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd);
static int asm_qud_rob(
    struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd);

static int slv(struct sim *sim, struct smtx *m, struct vec *b, struct vec *q)
{
    for (int i = 0; i < sim->msh->hxd.len; ++i)
        if (asm_hxd(sim, m, b, q, &sim->msh->hxd.dat[i]))
            return -1;

    struct ilog dir;

    if (log_new(&dir))
        return -1;

    for (int i = 0; i < sim->msh->qud.len; ++i) {
        struct qud     *qud = &sim->msh->qud.dat[i];
        struct bnd     *bnd = &sim->bnd.dat[qud->pid];
        struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

        switch (cnd->type) {
            case CND_BND_DIR:
                if (log_add(&dir, i)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_NEU:
                if (asm_qud_neu(sim, b, q, qud, cnd)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
            case CND_BND_ROB:
                if (asm_qud_rob(sim, m, b, q, qud, cnd)) {
                    log_cls(&dir);
                    return -1;
                }

                break;
        }
    }

    log_rst(&dir);

    for (int i = 0; !log_adv(&dir, &i);) {
        struct qud     *qud = &sim->msh->qud.dat[i];
        struct bnd     *bnd = &sim->bnd.dat[qud->pid];
        struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

        if (asm_qud_dir(sim, m, b, q, qud, cnd)) {
            log_cls(&dir);
            return -1;
        }
    }

    log_cls(&dir);

    if (errno != ENOENT)
        return -1;

    errno = 0;

    switch (sim->ops.ell.ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(m, q, b, sim->ops.ell.ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int asm_mov_mtx(struct smtx *m, int i, int j, double v);

static int asm_hxd(struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct hxd *hxd)
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

    struct obj *obj = &sim->obj.dat[hxd->pid];
    struct mat *mat = &sim->mat.dat[obj->mat];
    struct val *val = &sim->src.dat[obj->src];

    struct vtx *vtx = sim->msh->vtx.dat;

    if (mat->lam.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            lam[k] = mat->lam.as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

    if (mat->gam.type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            gam[k] = mat->gam.as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

    if (val->type == VAL_FUN)
        for (int k = 0; k < 8; ++k)
            src[k] = val->as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

    double hx = fabs(vtx[hxd->vtx[1]].x - vtx[hxd->vtx[0]].x);
    double hy = fabs(vtx[hxd->vtx[2]].y - vtx[hxd->vtx[0]].y);
    double hz = fabs(vtx[hxd->vtx[4]].z - vtx[hxd->vtx[0]].z);

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

        int mui = MU[DM[i]];
        int nui = NU[DM[i]];
        int tti = TT[DM[i]];

        for (int j = 0; j < 8; ++j) {
            int gj = hxd->vtx[j];

            int muj = MU[DM[j]];
            int nuj = NU[DM[j]];
            int ttj = TT[DM[j]];

            double mij = 0;

            if (mat->lam.type == VAL_FUN)
                for (int k = 0; k < 8; ++k) {
                    int muk = MU[DM[k]];
                    int nuk = NU[DM[k]];
                    int ttk = TT[DM[k]];

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
                    int muk = MU[DM[k]];
                    int nuk = NU[DM[k]];
                    int ttk = TT[DM[k]];

                    mij += gam[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                }
            else
                mij += mat->gam.as.num * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);

            asm_mov_mtx(m, gi, gj, mij);
        }

        double bi = 0;

        if (val->type == VAL_FUN)
            for (int k = 0; k < 8; ++k) {
                int muk = MU[DM[k]];
                int nuk = NU[DM[k]];
                int ttk = TT[DM[k]];

                bi += src[k] * (mx[muk][mui] * my[nuk][nui] * mz[ttk][tti]);
            }
        else
            bi = val->as.num * hx * hy * hz / 8;

        b->dat[gi] += bi;
    }

    return 0;
}

static int qud_nrm(struct qud *q, struct vtx *v)
{
    struct vtx *a = &v[q->vtx[0]];
    struct vtx *b = &v[q->vtx[1]];
    struct vtx *c = &v[q->vtx[2]];

    if (a->x == c->x && b->x == c->x)
        return 0;

    if (a->y == c->y && b->y == c->y)
        return 1;

    return 2;
}

static int asm_qud_dir(
    struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd)
{
    for (int i = 0; i < 4; ++i) {
        int gi = qud->vtx[i];

        m->dr[gi] = C;

        if (cnd->pps.dir.tgt.type == VAL_FUN)
            b->dat[gi] = C * cnd->pps.dir.tgt.as.fun(sim, gi, q->dat[gi]);
        else
            b->dat[gi] = C * cnd->pps.dir.tgt.as.num;
    }

    return 0;
}

static int asm_qud_neu(struct sim *sim, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd)
{
    static double tta[4];

    static double mxi[2][2];
    static double mzt[2][2];

    if (cnd->pps.neu.tta.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            tta[i] = cnd->pps.neu.tta.as.fun(sim, qud->vtx[i], q->dat[qud->vtx[i]]);

    int    nrm = qud_nrm(qud, sim->msh->vtx.dat);
    double hxi;
    double hzt;

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

        int mui = MU[DM[i]];
        int nui = NU[DM[i]];

        double bi = 0;

        if (cnd->pps.neu.tta.type == VAL_FUN)
            for (int k = 0; k < 4; ++k) {
                int muk = MU[DM[k]];
                int nuk = NU[DM[k]];

                bi += tta[k] * mxi[muk][mui] * mzt[nuk][nui];
            }
        else
            bi = cnd->pps.neu.tta.as.num * hxi * hzt / 4;

        b->dat[gi] += bi;
    }

    return 0;
}

static int asm_qud_rob(
    struct sim *sim, struct smtx *m, struct vec *b, struct vec *q, struct qud *qud, struct cnd_bnd *cnd)
{
    static double bet[4];
    static double ext[4];

    static double mxi[2][2];
    static double mzt[2][2];

    static double mnxi[2][2][2];
    static double mnzt[2][2][2];

    if (cnd->pps.rob.bet.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            bet[i] = cnd->pps.rob.bet.as.fun(sim, qud->vtx[i], q->dat[qud->vtx[i]]);

    if (cnd->pps.rob.ext.type == VAL_FUN)
        for (int i = 0; i < 4; ++i)
            ext[i] = cnd->pps.rob.ext.as.fun(sim, qud->vtx[i], q->dat[qud->vtx[i]]);

    int    nrm = qud_nrm(qud, sim->msh->vtx.dat);
    double hxi;
    double hzt;

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

            for (int k = 0; k < 2; ++k) {
                mnxi[i][j][k] = MN[i][j][k] * hxi;
                mnzt[i][j][k] = MN[i][j][k] * hzt;
            }
        }

    switch (cnd->pps.rob.bet.type) {
        case VAL_NUM:
            for (int i = 0; i < 4; ++i) {
                int gi = qud->vtx[i];
                int mui = MU[DM[i]];
                int nui = NU[DM[i]];

                for (int j = 0; j < 4; ++j) {
                    int gj = qud->vtx[j];
                    int muj = MU[DM[j]];
                    int nuj = NU[DM[j]];

                    double mij = cnd->pps.rob.bet.as.num * (mxi[muj][mui] * mzt[nuj][nui]);

                    if (asm_mov_mtx(m, gi, gj, mij))
                        return -1;
                }
            }

            switch (cnd->pps.rob.ext.type) {
                case VAL_NUM:
                    for (int i = 0; i < 4; ++i)
                        b->dat[qud->vtx[i]] += cnd->pps.rob.bet.as.num * cnd->pps.rob.ext.as.num * hxi * hzt / 4;

                    break;
                case VAL_FUN:
                    for (int i = 0; i < 4; ++i) {
                        int mui = MU[DM[i]];
                        int nui = NU[DM[i]];

                        double bi = 0;

                        for (int k = 0; k < 4; ++k) {
                            int muk = MU[DM[k]];
                            int nuk = NU[DM[k]];

                            bi += ext[k] * mxi[muk][mui] * mzt[nuk][nui];
                        }

                        b->dat[qud->vtx[i]] += cnd->pps.rob.bet.as.num * bi;
                    }

                    break;
            }

            break;
        case VAL_FUN:
            for (int i = 0; i < 4; ++i) {
                int gi = qud->vtx[i];
                int mui = MU[DM[i]];
                int nui = NU[DM[i]];

                for (int j = 0; j < 4; ++j) {
                    int gj = qud->vtx[j];
                    int muj = MU[DM[j]];
                    int nuj = NU[DM[j]];

                    double mij = 0;

                    for (int k = 0; k < 4; ++k) {
                        int muk = MU[DM[k]];
                        int nuk = NU[DM[k]];

                        mij += bet[k] * (mnxi[muk][muj][mui] * mnzt[nuk][nuj][nui]);
                    }

                    if (asm_mov_mtx(m, gi, gj, mij))
                        return -1;
                }
            }

            switch (cnd->pps.rob.ext.type) {
                case VAL_NUM:
                    for (int i = 0; i < 4; ++i) {
                        int mui = MU[DM[i]];
                        int nui = NU[DM[i]];

                        double bi = 0;

                        for (int k = 0; k < 4; ++k) {
                            int muk = MU[DM[k]];
                            int nuk = NU[DM[k]];

                            bi += bet[k] * mxi[muk][mui] * mzt[nuk][nui];
                        }

                        b->dat[qud->vtx[i]] += cnd->pps.rob.ext.as.num * bi;
                    }

                    break;
                case VAL_FUN:
                    for (int i = 0; i < 4; ++i) {
                        int mui = MU[DM[i]];
                        int nui = NU[DM[i]];

                        double bi = 0;

                        for (int k = 0; k < 4; ++k) {
                            int muk = MU[DM[k]];
                            int nuk = NU[DM[k]];

                            for (int j = 0; j < 4; ++j) {
                                int muj = MU[DM[j]];
                                int nuj = NU[DM[j]];

                                bi += bet[k] * ext[j] * mnxi[muk][muj][mui] * mnzt[nuk][nuj][nui];
                            }
                        }

                        b->dat[qud->vtx[i]] += bi;
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
