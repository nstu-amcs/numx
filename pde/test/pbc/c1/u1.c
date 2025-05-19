#include <numx/pde/sim.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5 * x * t * t * t - 2 * z * t * t + y * t + 1;
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 30 * x * t * t - 8 * z * t + 2 * y;
}

double initial(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->ops.tdd.beg;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5 * x * t * t * t - 2 * z * t * t + y * t + 1;
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -5 * t * t * t;
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return t;
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -2 * t * t;
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5 * x * t * t * t - 2 * z * t * t + y * t + 1 - 0.1 * t;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5 * x * t * t * t - 2 * z * t * t + y * t + 1 + 0.5 * t * t * t;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5 * x * t * t * t - 2 * z * t * t + y * t + 1 + 0.2 * t * t;
}
