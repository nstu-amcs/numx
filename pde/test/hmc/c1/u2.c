#include <numx/pde/sim.h>
#include <math.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * x * y * sin(f * t) + z * cos(f * t);
}

double ssrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -f * 2 * z - 2 * y;
}

double csrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * f * x * x * y;
}

double sneu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -2 * x * y;
}

double cneu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double sneu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * x;
}

double cneu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double sneu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double cneu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 1;
}

double srob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * x * y - 0.1 * x * x;
}

double crob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return z;
}

double srob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * x * y + 0.2 * x * y;
}

double crob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return z;
}

double srob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * x * y;
}

double crob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return z - 0.1;
}
