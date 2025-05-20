#include <numx/pde/sim.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z;
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0.4 * (x + y + z);
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return -5;
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5;
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 5;
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z + 0.5;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z + 0.5;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z - 0.5;
}

double rob4(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z - 0.5;
}

double rob5(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return x + y + z - 0.5;
}
