#include <numx/pde/sim.h>
#include <math.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return (x + y) * sin(f * t) + 2 * z * cos(f * t);
}

double ssrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -4 * z * f - f * f * (x + y);
}

double csrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * f * (x + y) - f * f * 2 * z;
}

double sneu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -1;
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

    return 1;
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

    return 2;
}

double srob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y - 0.1;
}

double crob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * z;
}

double srob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + 0.1;
}

double crob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * z;
}

double srob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y;
}

double crob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * z - 0.2;
}
