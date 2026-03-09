#include <stdlib.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>

#include "std.h"

static int fem_ctx_new(struct sim *sim, struct fem_std_ctx *ctx);
static int fem_ctx_cls(struct sim *sim, struct fem_std_ctx *ctx);
static int fem_slv(struct sim *sim, struct fem_std_ctx *ctx);

int fem_std_exe(struct sim *sim)
{
    struct fem_std_ctx ctx;

    if (fem_ctx_new(sim, &ctx)) {
        return -1;
    }

    if (fem_slv(sim, &ctx)) {
        fem_ctx_cls(sim, &ctx);
        return -1;
    }

    if (fem_ctx_cls(sim, &ctx)) {
        return -1;
    }

    return 0;
}

static inline void mtx_prep(struct smtx *mtx)
{
    memset(mtx, 0, sizeof(struct smtx));
}

static inline void vec_prep(struct vec *vec)
{
    memset(vec, 0, sizeof(struct vec));
}

static int fem_ctx_new(struct sim *sim, struct fem_std_ctx *ctx)
{
    mtx_prep(&ctx->mtx);
    mtx_prep(&ctx->sig);
    mtx_prep(&ctx->chi);
    vec_prep(&ctx->vec);
    vec_prep(&ctx->w0);
    vec_prep(&ctx->w1);
    vec_prep(&ctx->w2);
    vec_prep(&ctx->w3);

    int n = sim->msh->vtx.v2d.len;
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

    switch (sim->msh->type) {
        case MSH_C2D:
            for (int i = 0; i < sim->msh->qud.len; ++i) {
                int *vtx = sim->msh->qud.dat[i].vtx;

                for (int j = 0; j < 4; ++j)
                    for (int k = 0; k < 4; ++k)
                        if (vtx[k] < vtx[j] && !log_add(&map[vtx[j]], vtx[k]))
                            z += 1;
            }

            break;
        case MSH_C3D:
            for (int i = 0; i < sim->msh->hxd.len; ++i) {
                int *vtx = sim->msh->hxd.dat[i].vtx;

                for (int j = 0; j < 8; ++j)
                    for (int k = 0; k < 8; ++k)
                        if (vtx[k] < vtx[j] && !log_add(&map[vtx[j]], vtx[k]))
                            z += 1;
            }

            break;
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

            __attribute__((fallthrough));
        case SIM_PBC:
            if ((r = mtx_new(&ctx->sig, ctx->mtx.pps)))
                goto end;

            mtx_sdup(&ctx->mtx, &ctx->sig);

            __attribute__((fallthrough));
        case SIM_ELL:
            break;
    }

end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}

static int fem_ctx_cls(struct sim *sim, struct fem_std_ctx *ctx)
{
    mtx_cls(&ctx->mtx);
    vec_cls(&ctx->vec);

    switch (sim->eqn) {
        case SIM_HYP:
            mtx_cls(&ctx->chi);

            __attribute__((fallthrough));
        case SIM_PBC:
            mtx_cls(&ctx->sig);

            __attribute__((fallthrough));
        case SIM_ELL:
            break;
    }

    return 0;
}

static apx_fun _apx[] = {
    [FEM_BFS_LIN] = fem_std_lin_apx,
};

static fem_asm_fun _asm[] = {
    [FEM_BFS_LIN] = fem_std_lin_asm,
};

static int fem_ell_slv(struct sim *sim, struct fem_std_ctx *ctx);
static int fem_pbc_slv(struct sim *sim, struct fem_std_ctx *ctx);
static int fem_hyp_slv(struct sim *sim, struct fem_std_ctx *ctx);

static int fem_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    int r = 0;

    if ((r = sim->ops.exp.ini(sim))) {
        return r;
    }

    struct fem *fem = (struct fem *)sim->slv;

    sim->slv->apx = _apx[fem->ops.bfs];

    switch (sim->eqn) {
        case SIM_ELL:
            return fem_ell_slv(sim, ctx);
        case SIM_PBC:
            return fem_pbc_slv(sim, ctx);
        case SIM_HYP:
            return fem_hyp_slv(sim, ctx);
    }

    return -1;
}

static int fem_sys_slv(struct sim *sim, struct fem_std_ctx *ctx);

static int fem_ell_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    int r = 0;

    if ((r = vec_new(&ctx->w0, ctx->vec.n))) {
        goto fem_ell_slv_end;
    }

    sim->slv->run.bs = 1;
    sim->slv->run.wgt[0] = &ctx->w0;

    if ((r = fem_sys_slv(sim, ctx))) {
        goto fem_ell_slv_end;
    }

    if (sim->slv->itr_cbk.run) {
        sim->slv->itr_cbk.run(sim->slv->itr_cbk.ctx, sim);
    }

    if (sim->ops.exp.put_v) {
        sim->ops.exp.put_v(sim, sim->ops.exp.sol, sim->slv->run.wgt[0]);
    }

fem_ell_slv_end:
    vec_cls(&ctx->w0);
    return r;
}

static int fem_pbc_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    return -1;
}

static int fem_hyp_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    return -1;
}

static int fem_sys_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    struct fem     *fem = (struct fem *)sim->slv;
    struct slv_ops *ops = &sim->slv->ops;

    // if (ops->non.map)
    //     return slv_non(sim, ctx);

    if (_asm[fem->ops.bfs](sim, ctx)) {
        return -1;
    }

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