#include <numx/com/cmp.h>
#include <numx/com/log.h>

#include "lin/fem.h"

static int ctx_new(struct sim *sim, struct fem_std_ctx *ctx);
static int ctx_cls(struct sim *sim, struct fem_std_ctx *ctx);

int fem_std_exe(struct sim *sim)
{
    struct fem_std_ctx ctx;

    if (ctx_new(sim, &ctx))
        return -1;

    struct fem *fem = (struct fem *)sim->slv;

    switch (fem->ops.bss) {
        case FEM_BSS_LIN:
            return fem_std_lin_slv(sim, &ctx);
    }

    if (ctx_cls(sim, &ctx))
        return -1;

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

static int ctx_new(struct sim *sim, struct fem_std_ctx *ctx)
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

    switch (sim->eqn) {
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

static int ctx_cls(struct sim *sim, struct fem_std_ctx *ctx)
{
    mtx_cls(&ctx->mtx);
    vec_cls(&ctx->vec);

    switch (sim->eqn) {
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
