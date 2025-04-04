#include <numx/pde/sim.h>

double source(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;
    double z = v->z;

    return 0.4 * (5 + 0.2 * x + y + 30 * z + 0.5 * x * y + x * z + 10 * y * z + x * y * z);
}

double neu1(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return -1 - 2.5 * y - 5 * z - 5 * y * z;
}

double neu2(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 5 + 2.5 * x + 50 * z + 5 * x * z;
}

double neu3(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return 150 + 5 * x + 50 * y + 5 * x * y;
}

double rob1(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 9.1 + 11.25 * y + 50.5 * z + 30.5 * y * z;
}

double rob2(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 4.5 - 0.05 * x + 25 * z + 0.5 * x * z;
}

double rob3(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return -10 - 0.3 * x - 4 * y;
}
