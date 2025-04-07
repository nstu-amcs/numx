#include <numx/pde/sim.h>

double source1(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;
    double z = v->z;

    return 8 * x * (y * z + z + y + 1);
}

double source2(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    return 0;
}

double source3(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;
    double z = v->z;

    return 3 * x * z * (y + 1);
}

double neu1(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return -2 * x * (z + 1);
}

double neu2(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return -2 * x * z;
}

double neu3(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return -2 * (1 + x + z + x * z);
}

double neu4(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 2 * (1 + y + z + y * z);
}

double neu5(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return 2 * (1 + x + y + x * y);
}

double neu6(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double y = v->y;
    double z = v->z;

    return 2 * z * (1 + y);
}

double neu7(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return 2 * x * (1 + y);
}

double neu8(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return -2 * x * (1 + y);
}

double neu9(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double y = v->y;

    return -2 * (1 + x + y + x * y);
}

double rob1(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 8 * x * (z + 1);
}

double rob2(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 11.5 * x * z;
}

double rob3(struct sim *sim, struct vec *wgt, int vtx)
{
    (void)sim;
    (void)wgt;
    (void)vtx;

    struct vtx *v = &sim->msh->vtx.dat[vtx];

    double x = v->x;
    double z = v->z;

    return 4.5 * (1 + x + z + x * z);
}
