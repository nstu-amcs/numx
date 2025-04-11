#include <assert.h>

#include <numx/pde/sim.h>
#include <numx/vec/vec.h>

#include "fem_lin.h"

int fem_lin_pbc_mtx_asm(struct sim *sim, struct vec *wgt)
{
    assert(sim);
    assert(wgt);

    sim->rt.wgt = wgt;

    static double sig[8];

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

    struct fem *fem = (struct fem *)sim->slv;
    struct vtx *vtx = sim->msh->vtx.dat;

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
        struct hxd *hxd = &sim->msh->hxd.dat[hi];
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

                if (mat->sig.type == VAL_FUN)
                    for (int k = 0; k < 8; ++k) {
                        int muk = MU[k];
                        int nuk = NU[k];
                        int ttk = TT[k];

                        mij += sig[k] * (mnx[muk][muj][mui] * mny[nuk][nuj][nui] * mnz[ttk][ttj][tti]);
                    }
                else
                    mij = sig_n * (mx[muj][mui] * my[nuj][nui] * mz[ttj][tti]);

                asm_mov_mtx(&fem->prv.pbc, gi, gj, mij);
            }
        }
    }

    return 0;
}

static int slv_i2s(struct sim *sim);
static int slv_i3s(struct sim *sim);
static int slv_i4s(struct sim *sim);

int fem_lin_pbc_slv(struct sim *sim)
{
    assert(sim);

    switch (((struct fem *)sim->slv)->ops.tdd) {
        case FEM_TDD_I2S:
            return slv_i2s(sim);
        case FEM_TDD_I3S:
            return slv_i3s(sim);
        case FEM_TDD_I4S:
            return slv_i4s(sim);
    }

    return 0;
}

int fem_lin_pbc_i3s_asm(struct sim *sim, struct vec *w1, struct vec *w2, struct vec *w3);
int fem_lin_pbc_i4s_asm(struct sim *sim, struct vec *w1, struct vec *w2, struct vec *w3, struct vec *w4);

static int slv_i1s_uni(struct sim *sim, int t, struct vec *wgt);
static int slv_i2s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2);
static int slv_i3s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2, struct vec *w3);
static int slv_i4s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2, struct vec *w3, struct vec *w4);

static int slv_i2s(struct sim *sim)
{
    assert(sim);

    int r = 0;
    int n = sim->ops.tdd.num;
    int h = sim->ops.tdd.hop;

    assert(n > 0);

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  w1;
    struct vec  w2;

    if ((r = vec_new(&w1, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w2, fem->prv.vec.n)))
        goto end;

    if ((r = slv_i1s_uni(sim, 0, &w1)))
        goto end;

    for (int i = 0, t = h; i < n; i += 1, t += h) {
        if ((r = slv_i2s_uni(sim, t, &w1, &w2)))
            goto end;

        w1.dat = w2.dat;
    }

    return 0;

end:
    vec_cls(&w1);
    vec_cls(&w2);

    return r;
}

static int slv_i3s(struct sim *sim)
{
    assert(sim);

    int r = 0;
    int n = sim->ops.tdd.num;
    int h = sim->ops.tdd.hop;

    assert(n > 1);

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  w1;
    struct vec  w2;
    struct vec  w3;

    if ((r = vec_new(&w1, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w2, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w3, fem->prv.vec.n)))
        goto end;

    if ((r = slv_i1s_uni(sim, 0, &w1)))
        goto end;

    if ((r = slv_i2s_uni(sim, h, &w1, &w2)))
        goto end;

    for (int i = 1, t = 2 * h; i < n; i += 1, t += h) {
        if ((r = slv_i3s_uni(sim, t, &w1, &w2, &w3)))
            goto end;

        w1.dat = w2.dat;
    }

    return 0;

end:
    vec_cls(&w1);
    vec_cls(&w2);
    vec_cls(&w3);

    return r;
}

static int slv_i4s(struct sim *sim)
{
    assert(sim);

    int r = 0;
    int n = sim->ops.tdd.num;
    int h = sim->ops.tdd.hop;

    assert(n > 2);

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  w1;
    struct vec  w2;
    struct vec  w3;
    struct vec  w4;

    if ((r = vec_new(&w1, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w2, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w3, fem->prv.vec.n)))
        goto end;

    if ((r = vec_new(&w4, fem->prv.vec.n)))
        goto end;

    if ((r = slv_i1s_uni(sim, 0, &w1)))
        goto end;

    if ((r = slv_i2s_uni(sim, h, &w1, &w2)))
        goto end;

    if ((r = slv_i3s_uni(sim, 2 * h, &w1, &w2, &w3)))
        goto end;

    for (int i = 2, t = 3 * h; i < n; i += 1, t += h) {
        if ((r = slv_i4s_uni(sim, t, &w1, &w2, &w3, &w4)))
            goto end;

        w1.dat = w2.dat;
        w2.dat = w3.dat;
        w3.dat = w4.dat;
    }

    return 0;

end:
    vec_cls(&w1);
    vec_cls(&w2);
    vec_cls(&w3);
    vec_cls(&w4);

    return r;
}

static int slv_i1s_uni(struct sim *sim, int t, struct vec *wgt)
{
    for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
        struct hxd     *h = &sim->msh->hxd.dat[hi];
        struct obj     *o = &sim->obj.dat[h->pid];
        struct cnd_ini *i = &sim->cnd_ini.dat[o->ini];

        switch (i->tgt.type) {
            case VAL_NUM:
                for (int vi = 0; vi < 8; ++vi)
                    wgt->dat[h->vtx[vi]] = i->tgt.as.num;

                break;
            case VAL_FUN:
                for (int vi = 0; vi < 8; ++vi)
                    wgt->dat[h->vtx[vi]] = i->tgt.as.fun(sim);

                break;
        }
    }

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim);

    return 0;
}

static int slv_i2s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem_lin_ell_asm(sim, w2))
        return -1;

    if (fem_lin_pbc_i2s_asm(sim, w1, w2))
        return -1;

    switch (fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.mtx, w2, &fem->prv.vec, fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    if (fem->slv.itr.run)
        fem->slv.itr.run(fem->slv.itr.ctx, sim, t, w2);

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim, t, w2);

    return 0;
}

static int slv_i3s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2, struct vec *w3)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem_lin_ell_asm(sim, w3))
        return -1;

    if (fem_lin_pbc_i3s_asm(sim, w1, w2, w3))
        return -1;

    switch (fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.mtx, w3, &fem->prv.vec, fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    if (fem->slv.itr.run)
        fem->slv.itr.run(fem->slv.itr.ctx, sim, t, w3);

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim, t, w3);

    return 0;
}

static int slv_i4s_uni(struct sim *sim, int t, struct vec *w1, struct vec *w2, struct vec *w3, struct vec *w4)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem_lin_ell_asm(sim, w4))
        return -1;

    if (fem_lin_pbc_i4s_asm(sim, w1, w2, w3, w4))
        return -1;

    switch (fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.mtx, w4, &fem->prv.vec, fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    if (fem->slv.itr.run)
        fem->slv.itr.run(fem->slv.itr.ctx, sim, t, w4);

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim, t, w4);

    return 0;
}
