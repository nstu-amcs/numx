#include <numx/com/cmp.h>
#include <numx/com/log.h>

#include "lin/fem.h"

static int ctx_new(struct sim *sim, struct fem_hmc_ctx *ctx);
static int ctx_cls(struct sim *sim, struct fem_hmc_ctx *ctx);

int fem_hmc_exe(struct sim *sim)
{
    struct fem_hmc_ctx ctx;

    if (ctx_new(sim, &ctx))
        return -1;

    struct fem *fem = (struct fem *)sim->slv;

    switch (fem->ops.bss) {
        case FEM_BSS_LIN:
            return fem_hmc_lin_slv(sim, &ctx);
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

static int ctx_new(struct sim *sim, struct fem_hmc_ctx *ctx)
{
    mtx_prep(&ctx->mtx);
    vec_prep(&ctx->vec);
    vec_prep(&ctx->wgt);

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

    if ((r = mtx_new(&ctx->mtx, ((struct smtx_pps){2 * n, 4 * z + n}))))
        goto end;

    for (int i = 0, e = 0; i < n; ++i) {
        ctx->mtx.ia[2 * i] = e;
        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e += 2) {
            ctx->mtx.ja[e] = 2 * j;
            ctx->mtx.ja[e + 1] = 2 * j + 1;
        }

        ctx->mtx.ia[2 * i + 1] = e;
        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e += 2) {
            ctx->mtx.ja[e] = 2 * j;
            ctx->mtx.ja[e + 1] = 2 * j + 1;
        }

        ctx->mtx.ja[e++] = 2 * i;
    }

    ctx->mtx.ia[2 * n] = z * 4 + n;

    if ((r = vec_new(&ctx->vec, 2 * n)))
        goto end;

end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}

static int ctx_cls(struct sim *sim, struct fem_hmc_ctx *ctx)
{
    (void)sim;

    mtx_cls(&ctx->mtx);
    vec_cls(&ctx->vec);

    return 0;
}
