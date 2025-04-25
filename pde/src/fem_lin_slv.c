#include "fem_lin.h"

#include <numx/non/dif.h>

double fem_lin_apx(struct apx_fun_ctx *ctx, struct vec *vtx)
{
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

int fem_lin_slv(struct sim *sim)
{
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

static int ell_slv(struct sim *sim)
{
    int r = 0;

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  wgt;

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

static int pbc_slv(struct sim *sim)
{
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

    if (fem->slv.ops.tdd > 2) {
        if ((r = vec_new(&w2, fem->prv.vec.n)))
            goto end;

        fem->slv.run.tb = 3;
    }

    if (fem->slv.ops.tdd > 3) {
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
    int    num = sim->ops.tdd.num;

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

static int pbc_i1s_slv(struct sim *sim)
{
    struct fem *fem = (struct fem *)sim->slv;

    for (int i = 0; i < sim->msh->hxd.len; ++i) {
        struct hxd     *hxd = &sim->msh->hxd.dat[i];
        struct obj     *obj = &sim->obj.dat[hxd->pid];
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
                fem->slv.run.wgt[1]->dat[hxd->vtx[k]] = ini->tgt.as.fun(&ctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
            }
        } else {
            for (int k = 0; k < 8; ++k)
                fem->slv.run.wgt[1]->dat[hxd->vtx[k]] = ini->tgt.as.num;
        }
    }

    return 0;
}

static int slv_non(struct sim *sim);

static int sys_slv(struct sim *sim)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem->slv.ops.non.ops.fd)
        return slv_non(sim);

    if (fem_lin_asm(sim))
        return -1;

    switch (fem->slv.ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.mtx, fem->slv.run.wgt[0], &fem->prv.vec, &fem->slv.ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int slv_non(struct sim *sim)
{
    int r = 0;

    struct fem *fem = (struct fem *)sim->slv;
    struct vec *wgt = fem->slv.run.wgt[0];

    struct vec prv;
    struct vec acc;
    struct vec new;

    if ((r = vec_new(&acc, fem->prv.vec.n)))
        goto end;

    bool rlx = sim->slv->ops.non.ops.rlx;

    int    max = fem->slv.ops.non.ops.max;
    double err = fem->slv.ops.non.ops.err;
    double cur = 1;
    double nrm = 0;

    if ((r = asm_and_est(sim, &tmp, &cur)))
        goto end;

    vec_nrm(&fem->prv.vec, &nrm);

    for (int i = 0; i < max && cur > err; ++i) {
        if (i > 0) {
            switch (fem->slv.ops.iss.mod) {
                case ISS_BCG:
                    if ((r = iss_bcg_slv(&fem->prv.mtx, wgt, &fem->prv.vec, &fem->slv.ops.iss.ops.bcg)))
                        goto end;

                    break;
                default:
                    errno = ENOTSUP;
                    r = -1;
                    goto end;
            }

            smtx_rst(&fem->prv.mtx);
            vec_rst(&fem->prv.vec);

            if (rlx) {
                if ((r = rlx_and_est(sim, &prv, &tmp, &t2, &cur, 0.01, 0, 2)))
                    goto end;
            } else {
                if ((r = asm_and_est(sim, &tmp, &cur)))
                    goto end;
            }
        } else {
            if ((r = asm_and_est(sim, &tmp, &cur)))
                goto end;

            vec_nrm(&fem->prv.vec, &nrm);
        }

        cur = cur / nrm;

        fem->slv.ops.non.run.itr = i;
        fem->slv.ops.non.run.err = cur;

        if (fem->slv.ops.non.ops.itr.run)
            fem->slv.ops.non.ops.itr.run(fem->slv.ops.non.ops.itr.ctx, sim);

        if (rlx) {
            double *t = wgt->dat;
            wgt->dat = prv.dat;
            prv.dat = t;
        }
    }

end:
    vec_cls(&tmp);

    return r;
}

static int slv_new_asm(struct sim *sim);

static int slv_new(struct sim *sim)
{
    (void)sim;
    int r = 0;

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  tmp;

    if ((r = vec_new(&tmp, fem->prv.vec.n)))
        goto end;

    int    max = fem->slv.ops.non.ops.max;
    double err = fem->slv.ops.non.ops.err;
    double cur = 1;
    double tnm = 0;
    double bnm = 0;

    for (int i = 0; i < fem->prv.vec.n; ++i)
        fem->slv.run.wgt[0]->dat[i] = 1;

    for (int i = 0; i < max && cur > err; ++i) {
        if (i > 0) {
            if (slv_new_asm(sim))
                return -1;

            switch (fem->slv.ops.iss.mod) {
                case ISS_BCG:
                    if ((r = iss_bcg_slv(&fem->prv.mtx, fem->slv.run.wgt[0], &fem->prv.vec, &fem->slv.ops.iss.ops.bcg)))
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

        fem->slv.ops.non.run.itr = i;
        fem->slv.ops.non.run.err = cur;

        if (fem->slv.ops.non.ops.itr.run)
            fem->slv.ops.non.ops.itr.run(fem->slv.ops.non.ops.itr.ctx, sim);
    }

end:
    vec_cls(&tmp);

    return r;
}

static void dif_twk_wgt(void *ctx, double hop, struct dif_ops *ops)
{
    ((struct sim_fun_ctx *)ctx)->sim->slv->run.wgt[0]->dat[ops->var] += hop;
}

static int slv_new_asm(struct sim *sim)
{
    static double gnx[2][2];
    static double gny[2][2];
    static double gnz[2][2];

    static double mnx[2][2][2];
    static double mny[2][2][2];
    static double mnz[2][2][2];

    struct fem *fem = (struct fem *)sim->slv;
    struct vec *vtx = sim->msh->vtx.dat;

    struct sim_fun_ctx ctx = {.sim = sim};
    struct dif_ops     dif = {
            .hop = 0,
            .twk = dif_twk_wgt,
    };

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
                for (int i = 0; i < 8; ++i)
                    dlam[i] = mat->lam.ops.dif(&ctx, &vtx[hxd->vtx[i]]);

                break;
            case DIF_NUM:
                for (int i = 0; i < 8; ++i) {
                    int gi = hxd->vtx[i];

                    dif.var = gi;
                    dif.vtx = &vtx[gi];

                    dlam[i] = dif_tpm(&ctx, (double (*)(void *, struct vec *))mat->lam.as.fun, &dif);
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

                    mij += sim->slv->run.wgt[0]->dat[gk] *
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

                    bij += sim->slv->run.wgt[0]->dat[gk] *
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

                bi += bij * sim->slv->run.wgt[0]->dat[gj];

                if (smtx_add(&fem->prv.mtx, gi, gj, mij))
                    return -1;
            }

            fem->prv.vec.dat[gi] += bi;
        }
    }

    return 0;
}

static int rlx_and_est(
    struct sim *sim, struct vec *prv, struct vec *t1, struct vec *t2, double *err, double eps, double a, double b)
{
    double      del = eps / 2;
    struct vec *wgt = sim->slv->run.wgt[0];
    double     *dat = t1->dat;

    t1->dat = wgt->dat;
    wgt->dat = dat;

    while (b - a > eps) {
        double r1 = (a + b - del) / 2;
        double r2 = (a + b + del) / 2;
        double e1 = 0;
        double e2 = 0;

        vec_mul(t1, wgt, r1);
        vec_cmb(wgt, prv, wgt, 1 - r1);

        if (asm_and_est(sim, t2, &e1))
            return -1;

        vec_mul(t1, wgt, r2);
        vec_cmb(wgt, prv, wgt, 1 - r2);

        if (asm_and_est(sim, t2, &e2))
            return -1;

        if (e1 < e2)
            b = r2;
        else
            a = r1;
    }

    sim->slv->ops.non.run.rlx = (b - a) / 2;

    vec_mul(t1, wgt, sim->slv->ops.non.run.rlx);
    vec_cmb(wgt, prv, wgt, 1 - sim->slv->ops.non.run.rlx);

    if (asm_and_est(sim, t2, err))
        return -1;

    return 0;
}

static int asm_and_est(struct sim *sim, struct vec *tmp, double *err)
{
    if (fem_lin_asm(sim))
        return -1;

    struct fem *fem = (struct fem *)sim->slv;

    mtx_vmlt(&fem->prv.mtx, sim->slv->run.wgt[0], tmp);
    vec_cmb(tmp, &fem->prv.vec, tmp, -1);
    vec_nrm(tmp, err);

    return 0;
}
