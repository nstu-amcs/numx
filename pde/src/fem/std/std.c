#include <stdio.h>
#include <stdlib.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/non/opm.h>

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

    if ((r = smtx_new(&ctx->mtx, (struct smtx_pps){.n = n, .z = z})))
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

static apx_fun _dif[] = {
    [FEM_BFS_LIN] = fem_std_lin_dif,
};

static fem_asm_fun _asm[] = {
    [FEM_BFS_LIN] = fem_std_lin_asm,
};

static fem_asm_fun _new[] = {
    [FEM_BFS_LIN] = fem_std_lin_new,
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

    sim->slv->apx.run = _apx[fem->ops.bfs];
    sim->slv->apx.dif = _dif[fem->ops.bfs];

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

    if (sim->ops.exp.put_v) {
        sim->ops.exp.put_v(sim, sim->ops.exp.sol, 0, sim->slv->run.wgt[0]);
    }

    if (sim->slv->itr_cbk.run) {
        sim->slv->itr_cbk.run(sim->slv->itr_cbk.ctx, sim);
    }

fem_ell_slv_end:
    vec_cls(&ctx->w0);
    return r;
}

static int fem_pbc_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    (void)sim;
    (void)ctx;

    return -1;
}

static int fem_hyp_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    (void)sim;
    (void)ctx;

    return -1;
}

static int fem_sys_slv_non(struct sim *sim, struct fem_std_ctx *ctx);

static int fem_sys_slv(struct sim *sim, struct fem_std_ctx *ctx)
{
    int r = 0;

    struct fem     *fem = (struct fem *)sim->slv;
    struct slv_ops *ops = &sim->slv->ops;

    if (ops->non.enable) {
        return fem_sys_slv_non(sim, ctx);
    }

    if (_asm[fem->ops.bfs](sim, ctx)) {
        return -ENOTSUP;
    }

    printf("[fem][sys][asm] ok\n");
    vec_rst(&ctx->w0);

    switch (ops->iss.mod) {
        case ISS_BCG:
            if ((r = iss_bcg_slv(&ctx->mtx, &ctx->w0, &ctx->vec, &ops->iss.ops.bcg))) {
                return r;
            }

            break;
        case ISS_GMR:
            if ((r = iss_gmr_slv(&ctx->mtx, &ctx->w0, &ctx->vec, &ops->iss.ops.gmr))) {
                return r;
            }

            break;
        default:
            return -ENOTSUP;
    }

    printf("[fem][sys][slv] ok\n");

    return 0;
}

/**
 *  @brief Estimation context.
 */
struct est_ctx
{
    struct sim         *sim;
    struct fem_std_ctx *ctx; // Solver context.
    struct vec         *tmp; // Temporary vector to compute error.

    // Relaxation

    double      rlx; // Current relaxation.
    struct vec *prv;
    struct vec *upd;
};

/**
 *  @brief Estimate current system error.
 */
static double est(struct est_ctx *ctx, struct vec *wgt)
{
    struct vec est = {.n = wgt->n, .dat = wgt->dat};

    if (ctx->rlx != 0) {
        vec_mul(wgt, ctx->upd, ctx->rlx);
        vec_cmb(ctx->upd, ctx->prv, ctx->upd, 1 - ctx->rlx);

        est.dat = ctx->upd->dat;
    }

    // Assemble system with the current weights.

    vec_swp(&est, &ctx->ctx->w0);
    fem_std_lin_asm(ctx->sim, ctx->ctx);
    vec_swp(&est, &ctx->ctx->w0);

    // Estimate error of the current assembly.

    double err = 0;

    mtx_vmul(&ctx->ctx->mtx, &est, ctx->tmp);
    vec_cmb(ctx->tmp, &ctx->ctx->vec, ctx->tmp, -1);
    vec_nrm(ctx->tmp, &err);

    return err;
}

/**
 *  @brief Optimization tweak function to populate relaxation.
 */
static void est_twk_rlx(void *ctx, double val, struct opm_ops *ops)
{
    (void)ops;
    ((struct est_ctx *)ctx)->rlx = val;
}

/**
 *  @brief Nonlinear system solver.
 *
 *  Uses fixed-point iteration with Newton's linearization tweaks (if enabled).
 */
static int fem_sys_slv_non(struct sim *sim, struct fem_std_ctx *ctx)
{
    int r = 0;

    struct fem     *fem = (struct fem *)sim->slv;
    struct non_ops *ops = &sim->slv->ops.non.ops;

    struct vec tmp; // Temporary vector for error estimation.
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
        .ctx = ctx,
        .rlx = 0,
        .tmp = &tmp,
        .prv = &prv,
        .upd = &upd,
    };

    // Options for relaxation parameter optimization.
    struct opm_ops opm_ops = {
        .beg = 0.5,
        .end = 1.5,
        .eps = 0.01,
        .twk = est_twk_rlx,
        .var = -1,
        .vtx = &ctx->w0,
    };

    // Initialize field.

    if (ops->ini_cbk.run) {
        ops->ini_cbk.run(ops->ini_cbk.ctx, &ctx->w0);
    } else {
        vec_rst(&ctx->w0);
    }

    // Compute initial system error. Zero iteration.

    cur = est(&est_ctx, &ctx->w0);
    vec_nrm(&ctx->vec, &nrm);
    cur = cur / nrm;

    ops->run.err = cur;
    ops->run.itr = 0;
    ops->run.rlx = 0;

    if (ops->itr_cbk.run) {
        ops->itr_cbk.run(ops->itr_cbk.ctx, ops);
    }

    // Iteration process.

    for (int i = 1; i <= ops->max && cur > ops->err; ++i) {
        // If Newton's tweaks are enabled, assemble them.

        if (sim->slv->ops.non.ops.new) {
            if ((r = _new[fem->ops.bfs](sim, ctx))) {
                goto end;
            }
        }

        if (ops->rlx) {
            vec_swp(&ctx->w0, &prv);
        }

        vec_rst(&ctx->w0);

        switch (sim->slv->ops.iss.mod) {
            case ISS_BCG:
                if ((r = iss_bcg_slv(&ctx->mtx, &ctx->w0, &ctx->vec, &sim->slv->ops.iss.ops.bcg))) {
                    goto end;
                }

                break;
            case ISS_GMR:
                if ((r = iss_gmr_slv(&ctx->mtx, &ctx->w0, &ctx->vec, &sim->slv->ops.iss.ops.gmr))) {
                    goto end;
                }

                break;
            default:
                r = -ENOTSUP;
                goto end;
        }

        if (ops->rlx) {
            double opm = opm_loc_bis(&est_ctx, (double (*)(void *, struct vec *))est, &opm_ops);

            vec_mul(&ctx->w0, &upd, opm);
            vec_cmb(&upd, &prv, &upd, 1 - opm);
            vec_swp(&upd, &ctx->w0);

            est_ctx.rlx = 0;
            ops->run.rlx = opm;
        }

        // Compute error estimation for new solution.
        // Assebled system will be used in the next iteration.

        cur = est(&est_ctx, &ctx->w0) / nrm;

        ops->run.itr = i;
        ops->run.err = cur;

        if (ops->itr_cbk.run) {
            ops->itr_cbk.run(ops->itr_cbk.ctx, ops);
        }
    }

end:
    vec_cls(&tmp);

    if (ops->rlx) {
        vec_cls(&prv);
        vec_cls(&upd);
    }

    return r;
}
