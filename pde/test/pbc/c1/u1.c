#include <numx/pde/sim.h>
#include <math.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * y * sin(2 * t) + z * cos(2 * t);
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -4 * z * sin(2 * t) + 4 * x * y * cos(2 * t);
}

double initial(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    return target(ctx, vtx);
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -y * sin(2 * t);
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * sin(2 * t);
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return cos(2 * t);
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return (x * y - 0.1 * x) * sin(2 * t) + z * cos(2 * t);
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return (x * y + 0.1 * y) * sin(2 * t) + z * cos(2 * t);
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x * y * sin(2 * t) + (z - 0.1) * cos(2 * t);
}
