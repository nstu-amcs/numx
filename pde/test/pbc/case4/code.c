#include <numx/pde/sim.h>

double source(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;
    double z = v->z;

    int t = sim->rt.tv;

    return 0.4 * (x * y * z * t * t + 2 * t) + 2 * x * y * z * t + 2;
}

double initial(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    return 0;
}

double neu1(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    int t = sim->rt.tv;

    return -5 * y * z * t * t;
}

double neu2(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    int t = sim->rt.tv;

    return 5 * x * z * t * t;
}

double neu3(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    int t = sim->rt.tv;

    return 5 * x * y * t * t;
}

double rob1(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    int t = sim->rt.tv;

    return 20 * y * z * t * t + 2 * t + (y * z * t * t) / 2;
}

double rob2(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    int t = sim->rt.tv;

    return 2 * t - (x * z * t * t) / 2;
}

double rob3(struct sim *sim, int vtx)
{
    (void)sim;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    int t = sim->rt.tv;

    return 2 * t - (x * y * t * t) / 2;
}
