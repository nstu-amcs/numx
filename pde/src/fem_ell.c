#include <assert.h>
#include <stdlib.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>

static double fabs(double a)
{
    return a < 0 ? -a : a;
}

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

const int DM[8] = {2, 3, 7, 6, 0, 1, 5, 4};
const int RM[8] = {4, 5, 0, 1, 7, 6, 3, 2};

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

        for (int j = 0; log_adv(&map[i], &j); e++)
            m->ja[e] = j;

        log_cls(&map[i]);
    }

    m->ia[n] = z;

    free(map);

    return 0;
}

static int slv(struct sim *sim, struct smtx *m, struct vec *b, struct vec *q)
{
    double lam[8];
    double gam[8];
    double src[8];

    double gx[2][2];
    double gy[2][2];
    double gz[2][2];

    double mx[2][2];
    double my[2][2];
    double mz[2][2];

    double gnx[2][2];
    double gny[2][2];
    double gnz[2][2];

    double mnx[2][2][2];
    double mny[2][2][2];
    double mnz[2][2][2];

    struct vtx *vtx = sim->msh->vtx.dat;

    for (int h = 0; h < sim->msh->hxd.len; ++h) {
        struct hxd *hxd = &sim->msh->hxd.dat[h];

        struct obj *obj = &sim->obj.dat[hxd->pid];
        struct mat *mat = &sim->mat.dat[obj->mat];
        struct val *val = &sim->src.dat[obj->src];

        if (mat->lam.type == VAL_FUN)
            for (int k = 0; k < 8; ++k)
                lam[k] = mat->lam.as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

        if (mat->gam.type == VAL_FUN)
            for (int k = 0; k < 8; ++k)
                gam[k] = mat->gam.as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

        switch (val->type) {
            case VAL_NUM:
                for (int k = 0; k < 8; ++k)
                    src[k] = val->as.num;

                break;
            case VAL_FUN:
                for (int k = 0; k < 8; ++k)
                    src[k] = val->as.fun(sim, hxd->vtx[k], q->dat[hxd->vtx[k]]);

                break;
        }

        double hx = fabs(vtx[1].x - vtx[0].x);
        double hy = fabs(vtx[4].y - vtx[0].y);
        double hz = fabs(vtx[3].z - vtx[0].z);

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
            int mui = MU[DM[i]];
            int nui = NU[DM[i]];
            int tti = TT[DM[i]];

            for (int j = 0; j < 8; ++j) {
                int muj = MU[DM[j]];
                int nuj = NU[DM[j]];
                int ttj = TT[DM[j]];

                double tm = 0;
                double tb = 0;

                switch (mat->lam.type) {
                    case VAL_NUM:
                        tm = mat->lam.as.num *
                             (gx[mui][muj] * my[nui][nuj] * mz[tti][ttj] + mx[mui][muj] * gy[nui][nuj] * mz[tti][ttj] +
                                 mx[mui][muj] * my[nui][nuj] * gz[tti][ttj]);

                        break;
                    case VAL_FUN:
                        for (int k = 0; k < 8; ++k) {
                            int muk = MU[DM[k]];
                            int nuk = NU[DM[k]];
                            int ttk = TT[DM[k]];

                            tm += lam[k] * (gnx[mui][muj] * mny[nuk][nui][nuj] * mnz[ttk][tti][ttj] +
                                               mnx[muk][mui][muj] * gny[nui][nuj] * mnz[ttk][tti][ttj] +
                                               mnx[muk][mui][muj] * mny[nuk][nui][nuj] * gnz[tti][ttj]);
                        }

                        break;
                }

                switch (mat->gam.type) {
                    case VAL_NUM:
                        tm = mat->gam.as.num * mx[mui][muj] * my[nui][nuj] * mz[tti][ttj];
                        break;
                    case VAL_FUN:
                        for (int k = 0; k < 8; ++k) {
                            int muk = MU[DM[k]];
                            int nuk = NU[DM[k]];
                            int ttk = TT[DM[k]];

                            tm += gam[k] * mnx[muk][mui][muj] * mny[nuk][nui][nuj] * mnz[ttk][tti][ttj];
                        }

                        break;
                }

                tb = src[i] * mx[mui][muj] * my[nui][nuj] * mz[tti][ttj];

                int vi = hxd->vtx[i];
                int vj = hxd->vtx[j];

                b->dat[vi] += tb;

                if (vi > vj) {
                    int p = m->ia[vi];

                    while (m->ja[p] < vj)
                        ++p;

                    m->lr[p] += tm;

                    continue;
                }

                if (vi < vj) {
                    int p = m->ia[vj];

                    while (m->ja[p] < vi)
                        ++p;

                    m->ur[p] += tm;

                    continue;
                }

                m->dr[vi] += tm;
            }
        }
    }

    return 0;
}
