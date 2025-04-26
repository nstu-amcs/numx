#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>
#include <numx/pde/sim.h>

#include <prv/fem/lin.h>

int fem_new(struct fem *fem)
{
    assert(fem);

    fem->slv.exe = fem_exe;

    fem->ops.mod = FEM_STD;
    fem->ops.bss = FEM_BSS_LIN;

    fem->slv.ops.iss.mod = ISS_BCG;
    fem->slv.ops.iss.ops.bcg.con.sm = 0;
    fem->slv.ops.iss.ops.bcg.ops.err = 1e-10;
    fem->slv.ops.iss.ops.bcg.ops.itr.run = 0;
    fem->slv.ops.iss.ops.bcg.ops.max = 500;

    fem->slv.ops.non.mod = NON_FPI;
    fem->slv.ops.non.ops.itr.ctx = 0;
    fem->slv.ops.non.ops.itr.run = 0;
    fem->slv.ops.non.ops.max = 50;
    fem->slv.ops.non.ops.err = 1e-10;
    fem->slv.ops.non.ops.rlx = 1;

    fem->slv.ops.tdd = TDD_I2S;

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

static int fem_ctx_new(struct sim *sim, struct fem_ctx *ctx)
{
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
