#include <numx/pde/sim.h>
#include <math.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return pow(M_E, x) * sin(f * t);
}

double ssrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -pow(M_E, x);
}

double csrc(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 2 * f * pow(M_E, x);
}

double sneu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -pow(M_E, x);
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

    return 0;
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

    return 0;
}

double srob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return pow(M_E, x);
}

double crob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double srob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return pow(M_E, x) + 0.1 * pow(M_E, x);
}

double crob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double srob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return pow(M_E, x);
}

double crob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double f = ctx->sim->ops.hmc.frq;
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}
