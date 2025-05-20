#include <numx/pde/sim.h>

#define K 3

double target(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return x * sum;
}

double source(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = -1;
    double pow = 1;

    for (int p = 1; p <= K; ++p) {
        sum += sig * p * pow;
        sig *= -1;
        pow *= t;
    }

    return 2 * x * sum;
}

double initial(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return x * sum;
}

double neu1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return -sum;
}

double neu2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double neu3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    return 0;
}

double rob1(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return x * sum;
}

double rob2(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return x * sum + 0.1 * sum;
}

double rob3(struct sim_fun_ctx *ctx, struct vec *vtx)
{
    [[maybe_unused]] double t = ctx->sim->slv->run.tv;
    [[maybe_unused]] double x = vtx->dat[0];
    [[maybe_unused]] double y = vtx->dat[1];
    [[maybe_unused]] double z = vtx->dat[2];

    double sum = 0;
    double sig = 1;
    double pow = 1;

    for (int p = 0; p <= K; ++p) {
        sum += sig * pow;
        sig *= -1;
        pow *= t;
    }

    return x * sum;
}
