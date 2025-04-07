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

    return 0.4 * (x + y + z);
}

double neu1(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return -5;
}

double neu2(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 5;
}

double neu3(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 5;
}

double rob1(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 15.5 + y + z;
}

double rob2(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 20.5 + y + z;
}

double rob3(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return -0.5 + x + z;
}

double rob4(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 3.5 + x + z;
}

double rob5(struct sim *sim, struct vec* wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return -0.5 + x + y;
}
