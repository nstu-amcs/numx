#include <numx/pde/sim.h>

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 5 + 0.2 * x + y + 30 * z + 0.5 * x * y + x * z + 10 * y * z + x * y * z;
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 0.4 * (5 + 0.2 * x + y + 30 * z + 0.5 * x * y + x * z + 10 * y * z + x * y * z);
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return -1 - 2.5 * y - 5 * z - 5 * y * z;
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double z = vtx->dat[2];

    return 5 + 2.5 * x + 50 * z + 5 * x * z;
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];

    return 150 + 5 * x + 50 * y + 5 * x * y;
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return 9.1 + 11.25 * y + 50.5 * z + 30.5 * y * z;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double z = vtx->dat[2];

    return 4.5 - 0.05 * x + 25 * z + 0.5 * x * z;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];

    return -10 - 0.3 * x - 4 * y;
}
