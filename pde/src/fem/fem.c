#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>
#include <numx/pde/sim.h>

#include <prv/fem/lin.h>

int fem_new(struct fem *fem)
{
    assert(fem);

    fem->slv.ops.tdd = TDD_I2S;

    fem->slv.ops.non.mod = NON_FPI;
    fem->slv.ops.non.map = 0;
    fem->slv.ops.non.ops.max = 200;
    fem->slv.ops.non.ops.err = 1e-10;
    fem->slv.ops.non.ops.rlx = false;
    fem->slv.ops.non.ops.itr.ctx = NULL;
    fem->slv.ops.non.ops.itr.run = NULL;
    fem->slv.ops.non.ops.dif = DIF_NUM;
    fem->slv.ops.non.ops.run.itr = 0;
    fem->slv.ops.non.ops.run.err = 0;
    fem->slv.ops.non.ops.run.rlx = 0;

    fem->slv.ops.iss.mod = ISS_BCG;
    fem->slv.ops.iss.ops.bcg.ops.max = 500;
    fem->slv.ops.iss.ops.bcg.ops.err = 1e-10;
    fem->slv.ops.iss.ops.bcg.ops.itr.ctx = NULL;
    fem->slv.ops.iss.ops.bcg.ops.itr.run = NULL;
    fem->slv.ops.iss.ops.bcg.ops.run.itr = 0;
    fem->slv.ops.iss.ops.bcg.ops.run.err = 0;
    fem->slv.ops.iss.ops.bcg.con.sm = NULL;

    fem->slv.exe = fem_exe;
    fem->slv.apx = NULL;

    fem->slv.itr.ctx = NULL;
    fem->slv.itr.run = NULL;

    fem->slv.run.wgt[0] = NULL;
    fem->slv.run.wgt[1] = NULL;
    fem->slv.run.wgt[2] = NULL;
    fem->slv.run.wgt[3] = NULL;

    fem->slv.run.bs = 1;
    fem->slv.run.ti = 0;
    fem->slv.run.tv = 0;

    fem->ops.mod = FEM_STD;
    fem->ops.bss = FEM_BSS_LIN;

    return 0;
}

static int fem_ctx_new(struct sim *sim, struct fem_ctx *ctx);
static int fem_ctx_cls(struct sim *sim, struct fem_ctx *ctx);

int fem_exe(struct sim *sim)
{
    assert(sim);

    struct fem_ctx ctx;

    if (fem_ctx_new(sim, &ctx))
        return -1;

    if (sim->ops.exp.ini(sim))
        return -1;

    switch (((struct fem *)sim->slv)->ops.bss) {
        case FEM_BSS_LIN:
            return fem_lin_slv(sim, &ctx);
    }

    fem_ctx_cls(sim, &ctx);

    return 0;
}

static void mtx_prep(struct smtx *mtx)
{
    mtx->dr = NULL;
    mtx->lr = NULL;
    mtx->ur = NULL;
    mtx->ia = NULL;
    mtx->ja = NULL;
    mtx->pps.n = 0;
    mtx->pps.z = 0;
}

static void vec_prep(struct vec *vec)
{
    vec->n = 0;
    vec->dat = NULL;
}

static int fem_ctx_new(struct sim *sim, struct fem_ctx *ctx)
{
    mtx_prep(&ctx->mtx);
    mtx_prep(&ctx->sig);
    mtx_prep(&ctx->chi);

    vec_prep(&ctx->vec);
    vec_prep(&ctx->w0);
    vec_prep(&ctx->w1);
    vec_prep(&ctx->w2);
    vec_prep(&ctx->w3);

    int n = sim->msh->vtx.len;
    int z = 0;
    int r = 0;

    struct ilog *map = malloc(sizeof(struct ilog) * n);

    if (!map) {
        r = -1;
        goto end;
    }

    for (int i = 0; i < n; ++i) {
        if ((r = log_new(&map[i])))
            goto end;

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

    if ((r = mtx_new(&ctx->mtx, ((struct smtx_pps){n, z}))))
        goto end;

    for (int i = 0, e = 0; i < n; ++i) {
        ctx->mtx.ia[i] = e;

        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e++)
            ctx->mtx.ja[e] = j;
    }

    ctx->mtx.ia[n] = z;

    if ((r = vec_new(&ctx->vec, n)))
        goto end;

    switch (sim->mod) {
        case SIM_HYP:
            if ((r = mtx_new(&ctx->chi, ctx->mtx.pps)))
                goto end;

            mtx_sdup(&ctx->mtx, &ctx->chi);

            [[fallthrough]];
        case SIM_PBC:
            if ((r = mtx_new(&ctx->sig, ctx->mtx.pps)))
                goto end;

            mtx_sdup(&ctx->mtx, &ctx->sig);

            [[fallthrough]];
        case SIM_ELL:
            break;
    }

end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}

static int fem_ctx_cls(struct sim *sim, struct fem_ctx *ctx)
{
    mtx_cls(&ctx->mtx);
    vec_cls(&ctx->vec);

    switch (sim->mod) {
        case SIM_HYP:
            mtx_cls(&ctx->chi);

            [[fallthrough]];
        case SIM_PBC:
            mtx_cls(&ctx->sig);

            [[fallthrough]];
        case SIM_ELL:
            break;
    }

    return 0;
}
