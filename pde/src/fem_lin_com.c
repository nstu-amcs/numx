#include "fem_lin.h"

const double C = 10e10;

const double G[2][2] = {
    {1.0,  -1.0},
    {-1.0, 1.0 }
};

const double M[2][2] = {
    {2.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 6.0}
};

const double GN[2][2] = {
    {1.0 / 2.0,  -1.0 / 2.0},
    {-1.0 / 2.0, 1.0 / 2.0 }
};

const double MN[2][2][2] = {
    {{1.0 / 4.0, 1.0 / 12.0},  {1.0 / 12.0, 1.0 / 12.0}},
    {{1.0 / 12.0, 1.0 / 12.0}, {1.0 / 12.0, 1.0 / 4.0} }
};

const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

double fem_lin_apx(void *ctx, struct sim *sim, struct vec *wgt, int hxd, struct vtx *vtx)
{
    (void)ctx;

    struct vtx *v = sim->msh->vtx.dat;
    struct hxd *h = &sim->msh->hxd.dat[hxd];

    double *w = wgt->dat;

    int v0 = h->vtx[0];
    int v7 = h->vtx[7];

    double x1 = v[v0].x;
    double x2 = v[v7].x;
    double y1 = v[v0].y;
    double y2 = v[v7].y;
    double z1 = v[v0].z;
    double z2 = v[v7].z;

    double hm = (x2 - x1) * (y2 - y1) * (z2 - z1);

    double x = vtx->x;
    double y = vtx->y;
    double z = vtx->z;
    double r = 0;

    r += w[h->vtx[0]] * (x2 - x) * (y2 - y) * (z2 - z) / hm;
    r += w[h->vtx[1]] * (x - x1) * (y2 - y) * (z2 - z) / hm;
    r += w[h->vtx[2]] * (x2 - x) * (y - y1) * (z2 - z) / hm;
    r += w[h->vtx[3]] * (x - x1) * (y - y1) * (z2 - z) / hm;
    r += w[h->vtx[4]] * (x2 - x) * (y2 - y) * (z - z1) / hm;
    r += w[h->vtx[5]] * (x - x1) * (y2 - y) * (z - z1) / hm;
    r += w[h->vtx[6]] * (x2 - x) * (y - y1) * (z - z1) / hm;
    r += w[h->vtx[7]] * (x - x1) * (y - y1) * (z - z1) / hm;

    return r;
}

