#include "fem.h"

int fem_hmc_lin_slv(struct sim *sim, struct fem_hmc_ctx *ctx)
{
    if (vec_new(&ctx->wgt, ctx->vec.n))
        return -1;

    sim->slv->run.wgt[0] = &ctx->wgt;
    sim->slv->run.bs = 1;

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

    if (sim->slv->itr.run)
        sim->slv->itr.run(sim->slv->itr.ctx, sim);

    if (sim->ops.exp.put)
        sim->ops.exp.put(sim);

    return 0;
}
