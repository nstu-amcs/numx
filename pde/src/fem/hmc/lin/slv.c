#include "fem.h"
#include <math.h>

double fem_hmc_lin_apx(struct apx_fun_ctx *ctx, struct vec *vtx)
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

    double b1 = (x2 - x) * (y2 - y) * (z2 - z) / hm;
    double b2 = (x - x1) * (y2 - y) * (z2 - z) / hm;
    double b3 = (x2 - x) * (y - y1) * (z2 - z) / hm;
    double b4 = (x - x1) * (y - y1) * (z2 - z) / hm;
    double b5 = (x2 - x) * (y2 - y) * (z - z1) / hm;
    double b6 = (x - x1) * (y2 - y) * (z - z1) / hm;
    double b7 = (x2 - x) * (y - y1) * (z - z1) / hm;
    double b8 = (x - x1) * (y - y1) * (z - z1) / hm;

    double us = (w[h->vtx[0] * 2] * b1) +
                (w[h->vtx[1] * 2] * b2) +
                (w[h->vtx[2] * 2] * b3) +
                (w[h->vtx[3] * 2] * b4) +
                (w[h->vtx[4] * 2] * b5) +
                (w[h->vtx[5] * 2] * b6) +
                (w[h->vtx[6] * 2] * b7) +
                (w[h->vtx[7] * 2] * b8);

    double uc = (w[h->vtx[0] * 2 + 1] * b1) +
                (w[h->vtx[1] * 2 + 1] * b2) +
                (w[h->vtx[2] * 2 + 1] * b3) +
                (w[h->vtx[3] * 2 + 1] * b4) +
                (w[h->vtx[4] * 2 + 1] * b5) +
                (w[h->vtx[5] * 2 + 1] * b6) +
                (w[h->vtx[6] * 2 + 1] * b7) +
                (w[h->vtx[7] * 2 + 1] * b8);

    double t = ctx->sim->slv->run.tv;
    double f = ctx->sim->ops.hmc.frq;

    return us * sin(f * t) + uc * cos(f * t);
}

int fem_hmc_lin_slv(struct sim *sim, struct fem_hmc_ctx *ctx)
{
    sim->slv->apx = fem_hmc_lin_apx;

    if (vec_new(&ctx->wgt, ctx->vec.n))
        return -1;

    sim->slv->run.wgt[0] = &ctx->wgt;
    sim->slv->run.bs = 1;
    sim->slv->run.tv = sim->ops.tdd.beg;
    sim->slv->run.ti = 1;

    if (fem_hmc_lin_asm(sim, ctx))
        return -1;

    switch (sim->slv->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&ctx->mtx, &ctx->wgt, &ctx->vec,
                    &sim->slv->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    double hop = sim->ops.tdd.hop;
    int    num = sim->ops.tdd.num;

    if (sim->eqn == SIM_ELL)
        num = 1;

    for (int i = 1; i <= num; ++i) {
        if (sim->slv->itr.run)
            sim->slv->itr.run(sim->slv->itr.ctx, sim);

        if (sim->ops.exp.put)
            sim->ops.exp.put(sim);

        sim->slv->run.tv += hop;
        sim->slv->run.ti += 1;
    }

    vec_cls(&ctx->wgt);

    return 0;
}
