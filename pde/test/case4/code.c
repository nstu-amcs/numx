#include <numx/pde/sim.h>

double lam(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 10 * wgt->dat[vtx];
}

double source(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;
    double z = v->z;

    return -30 + 0.4 * (x + y + z);
}

double neu1(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return -10 * wgt->dat[vtx];
}

double neu2(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 10 * wgt->dat[vtx];
}

double neu3(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 10 * wgt->dat[vtx];
}

double rob1(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 20 + y + z + wgt->dat[vtx];
}

double rob2(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return x + z - wgt->dat[vtx];
}

double rob3(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return x + y - wgt->dat[vtx];
}
