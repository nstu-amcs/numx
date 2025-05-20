#include <numx/non/dif.h>
#include <numx/non/opm.h>

#include <prv/fem/lin.h>

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

static int ell_slv(struct sim *sim, struct fem_ctx *ctx);
static int pbc_slv(struct sim *sim, struct fem_ctx *ctx);
static int hyp_slv(struct sim *sim, struct fem_ctx *ctx);

int fem_lin_slv(struct sim *sim, struct fem_ctx *ctx)
{
    assert(sim);

    sim->slv->apx = fem_lin_apx;

    switch (sim->mod) {
        case SIM_ELL:
            return ell_slv(sim, ctx);
        case SIM_PBC:
            return pbc_slv(sim, ctx);
        case SIM_HYP:
            return hyp_slv(sim, ctx);
    }

    return 0;
}

static int sys_slv(struct sim *sim, struct fem_ctx *ctx);

static int ell_slv(struct sim *sim, struct fem_ctx *ctx)
{
    int r = 0;

    if ((r = vec_new(&ctx->w0, ctx->vec.n)))
        goto end;

    sim->slv->run.bs = 1;
    sim->slv->run.wgt[0] = &ctx->w0;

    if ((r = sys_slv(sim, ctx)))
        goto end;

    if (sim->slv->itr.run)
        sim->slv->itr.run(sim->slv->itr.ctx, sim);

    if (sim->ops.exp.put)
        sim->ops.exp.put(sim);

end:
    vec_cls(&ctx->w0);
    return r;
}

static int pbc_i1s_slv(struct sim *sim, struct fem_ctx *ctx);
static int pbc_ctx_shr(struct sim *sim, struct fem_ctx *ctx);

static int pbc_slv(struct sim *sim, struct fem_ctx *ctx)
{
    int r = 0;

    if ((r = vec_new(&ctx->w0, ctx->vec.n)))
        goto end;

    if ((r = vec_new(&ctx->w1, ctx->vec.n)))
        goto end;

    if ((r = vec_new(&ctx->tmp, ctx->vec.n)))
        goto end;

    sim->slv->run.bs = 2;

    if (sim->slv->ops.tdd > 2) {
        if ((r = vec_new(&ctx->w2, ctx->vec.n)))
            goto end;

        sim->slv->run.bs = 3;
    }

    if (sim->slv->ops.tdd > 3) {
        if ((r = vec_new(&ctx->w3, ctx->vec.n)))
            goto end;

        sim->slv->run.bs = 4;
    }

    sim->slv->run.wgt[0] = &ctx->w0;
    sim->slv->run.wgt[1] = &ctx->w1;
    sim->slv->run.wgt[2] = &ctx->w2;
    sim->slv->run.wgt[3] = &ctx->w3;

    double beg = sim->ops.tdd.beg;
    double hop = sim->ops.tdd.hop;
    int    num = sim->ops.tdd.num;

    sim->slv->run.tv = beg;

    for (int i = 0; i <= num; ++i) {
        sim->slv->run.ti = i;

        if (i < sim->slv->ops.ini.num) {
            pbc_i1s_slv(sim, ctx);
        } else {
            sys_slv(sim, ctx);
        }

        if (sim->slv->itr.run)
            sim->slv->itr.run(sim->slv->itr.ctx, sim);

        if (sim->ops.exp.put)
            sim->ops.exp.put(sim);

        sim->slv->run.tv += hop;

        pbc_ctx_shr(sim, ctx);
    }

end:
    vec_cls(&ctx->w0);
    vec_cls(&ctx->w1);
    vec_cls(&ctx->w2);
    vec_cls(&ctx->w3);

    vec_cls(&ctx->tmp);

    return r;
}

static int pbc_i1s_slv(struct sim *sim, struct fem_ctx *ctx)
{
    for (int i = 0; i < sim->msh->hxd.len; ++i) {
        struct hxd     *hxd = &sim->msh->hxd.dat[i];
        struct obj     *obj = &sim->obj.dat[hxd->pid];
        struct cnd_ini *ini = &sim->cnd_ini.dat[obj->ini];

        if (ini->tgt.type == VAL_FUN) {
            struct sim_fun_ctx fctx = {
                .sim = sim,
                .vtx = -1,
                .qud = -1,
                .hxd = i,
            };

            for (int k = 0; k < 8; ++k) {
                fctx.vtx = hxd->vtx[k];
                ctx->w0.dat[hxd->vtx[k]] =
                    ini->tgt.as.fun(&fctx, &sim->msh->vtx.dat[hxd->vtx[k]]);
            }
        } else {
            for (int k = 0; k < 8; ++k)
                ctx->w0.dat[hxd->vtx[k]] = ini->tgt.as.num;
        }
    }

    return 0;
}

static int pbc_ctx_shr(struct sim *sim, struct fem_ctx *ctx)
{
    struct vec *w0 = &ctx->w0;
    struct vec *w1 = &ctx->w1;
    struct vec *w2 = &ctx->w2;
    struct vec *w3 = &ctx->w3;

    switch (sim->slv->ops.tdd) {
        case TDD_I2S:
            vec_swp(w0, w1);
            break;
        case TDD_I3S:
            vec_swp(w1, w2);
            vec_swp(w0, w1);
            break;
        case TDD_I4S:
            vec_swp(w2, w3);
            vec_swp(w1, w2);
            vec_swp(w0, w1);
            break;
    }

    return 0;
}

static int slv_non(struct sim *sim, struct fem_ctx *ctx);

static int sys_slv(struct sim *sim, struct fem_ctx *ctx)
{
    struct slv_ops *ops = &sim->slv->ops;

    if (ops->non.map)
        return slv_non(sim, ctx);

    if (fem_lin_asm(sim, ctx))
        return -1;

    vec_rst(&ctx->w0);

    switch (ops->iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&ctx->mtx, &ctx->w0, &ctx->vec, &ops->iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

struct est_ctx
{
    struct sim *sim;
    struct vec *tmp;

    double      rlx;
    struct vec *prv;
    struct vec *upd;

    struct fem_ctx *ctx;
};

static double est(struct est_ctx *ctx, struct vec *wgt)
{
    struct vec est = {.n = wgt->n, .dat = wgt->dat};

    if (ctx->rlx != 0) {
        vec_mul(wgt, ctx->upd, ctx->rlx);
        vec_cmb(ctx->upd, ctx->prv, ctx->upd, 1 - ctx->rlx);

        est.dat = ctx->upd->dat;
    }

    vec_swp(&est, &ctx->ctx->w0);
    fem_lin_asm(ctx->sim, ctx->ctx);
    vec_swp(&est, &ctx->ctx->w0);

    double err = 0;

    mtx_vmul(&ctx->ctx->mtx, &est, ctx->tmp);
    vec_cmb(ctx->tmp, &ctx->ctx->vec, ctx->tmp, -1);
    vec_nrm(ctx->tmp, &err);

    return err;
}

void est_twk_rlx(void *ctx, double val, struct opm_ops *)
{
    ((struct est_ctx *)ctx)->rlx = val;
}

static int slv_non(struct sim *sim, struct fem_ctx *ctx)
{
    int r = 0;

    struct non_ops *ops = &sim->slv->ops.non.ops;

    struct vec tmp;
    struct vec prv;
    struct vec upd;

    double nrm = 0;
    double cur = 0;

    if ((r = vec_new(&tmp, ctx->vec.n)))
        goto end;

    if (ops->rlx) {
        if ((r = vec_new(&prv, ctx->vec.n)))
            goto end;

        if ((r = vec_new(&upd, ctx->vec.n)))
            goto end;
    }

    struct est_ctx est_ctx = {
        .sim = sim,
        .rlx = 0,
        .ctx = ctx,
        .tmp = &tmp,
        .prv = &prv,
        .upd = &upd,
    };

    struct opm_ops opm_ops = {
        .beg = 0.5,
        .end = 1.5,
        .eps = 0.01,
        .twk = est_twk_rlx,
        .var = -1,
        .vtx = &ctx->w0,
    };

    for (int i = 0; i < sim->msh->vtx.len; ++i) {
        struct vec *v = &sim->msh->vtx.dat[i];

        double x = v->dat[0];
        double y = v->dat[1];
        double z = v->dat[2];

        ctx->w0.dat[i] = 0.6 * x + 0.7 * y + 0.8 * z;
    }

    ops->run.itr = 0;
    ops->run.rlx = 0;
    cur = est(&est_ctx, &ctx->w0);

    vec_nrm(&ctx->vec, &nrm);

    cur = cur / nrm;
    ops->run.err = cur;

    if (ops->itr.run)
        ops->itr.run(ops->itr.ctx, ops);

    for (int i = 1; i <= ops->max && cur > ops->err; ++i) {
        ops->run.itr = i;

        if (sim->slv->ops.non.mod == NON_NEW)
            if (fem_lin_new(sim, ctx))
                return -1;

        if (ops->rlx)
            vec_swp(&ctx->w0, &prv);

        vec_rst(&ctx->w0);

        switch (sim->slv->ops.iss.mod) {
            case ISS_BCG:
                if ((r = iss_bcg_slv(&ctx->mtx, &ctx->w0, &ctx->vec,
                         &sim->slv->ops.iss.ops.bcg)))
                    goto end;

                break;
            default:
                errno = ENOTSUP;
                r = -1;
                goto end;
        }

        if (ops->rlx) {
            double opm = opm_loc_bis(
                &est_ctx, (double (*)(void *, struct vec *))est, &opm_ops);

            vec_mul(&ctx->w0, &upd, opm);
            vec_cmb(&upd, &prv, &upd, 1 - opm);
            vec_swp(&upd, &ctx->w0);

            est_ctx.rlx = 0;
            ops->run.rlx = opm;
        }

        cur = est(&est_ctx, &ctx->w0) / nrm;
        ops->run.err = cur;

        if (ops->itr.run)
            ops->itr.run(ops->itr.ctx, ops);
    }

end:
    vec_cls(&tmp);

    if (ops->rlx) {
        vec_cls(&prv);
        vec_cls(&upd);
    }

    return r;
}

static int hyp_slv(struct sim *, struct fem_ctx *)
{
    return 0;
}
