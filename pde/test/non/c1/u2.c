#include <numx/pde/sim.h>

#define U(x, y, z)  ((x) + (y) + (z))
#define L(x, y, z)  (U(x, y, z))
#define LX(x, y, z) (1)
#define LY(x, y, z) (1)
#define LZ(x, y, z) (1)
#define BET         10.0
#define GAM         0.4

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return U(x, y, z);
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return -LX(x, y, z) - LY(x, y, z) - LZ(x, y, z) + GAM * U(x, y, z);
}

double lam(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return x + y + z;
}

double dlam(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    return 0;
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return -L(x, y, z);
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return L(x, y, z);
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return L(x, y, z);
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return U(x, y, z) - L(x, y, z) / BET;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return U(x, y, z) + L(x, y, z) / BET;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    (void)ctx;
    (void)vtx;

    double x = vtx->dat[0];
    double y = vtx->dat[1];
    double z = vtx->dat[2];

    return U(x, y, z) - L(x, y, z) / BET;
}
