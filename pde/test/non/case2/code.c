#include <numx/pde/sim.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return x + y + z;
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return -30 + 0.4 * (x + y + z);
}

double lam(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    return 10 * ctx->sim->slv->run.wgt[0]->dat[ctx->vtx] + 1;
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return -10 * (x + y + z) - 1;
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 10 * (x + y + z) + 1;
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 10 * (x + y + z) + 1;
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 20.1 + 2 * y + 2 * z + x;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double y = vtx->dat[1];

    return -y - 0.1;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double z = vtx->dat[2];

    return -z - 0.1;
}
