#include <assert.h>
#include </usr/include/cgnslib.h>
#include <stdio.h>
#include <stdlib.h>

#include <numx/msh/umsh.h>

static int umsh_exp_gns_c2d(struct umsh *msh, const char* fname);
static int umsh_exp_gns_c3d(struct umsh *msh, const char* fname);

int umsh_exp_gns(struct umsh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    char fname[256];
    sprintf(fname, "%s/%s.cgns", dir, pfx);

    switch (msh->type) {
        case MSH_C2D:
            return umsh_exp_gns_c2d(msh, fname);
        case MSH_C3D:
            return umsh_exp_gns_c3d(msh, fname);
        default:
            return -1;
    }
}

static int umsh_exp_gns_c2d(struct umsh *msh, const char *fname)
{
    int fi;
    int bi;
    int zi;
    int ci;
    int si;

    int dim[3] = {msh->vtx.v2d.len, msh->qud.len, 0};

    double *x = malloc(sizeof(double) * msh->vtx.v2d.len);
    double *y = malloc(sizeof(double) * msh->vtx.v2d.len);

    for (int i = 0; i < msh->vtx.v2d.len; ++i) {
        x[i] = msh->vtx.v2d.dat[i].dat[0];
        y[i] = msh->vtx.v2d.dat[i].dat[1];
    }

    cg_open(fname, CG_MODE_WRITE, &fi);

    cg_base_write(fi, "Base", 2, 2, &bi);
    cg_zone_write(fi, bi, "Zone 1", dim, Unstructured, &zi);

    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateX", x, &ci);
    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateY", y, &ci);

    int *node = malloc(sizeof(int *) * msh->qud.len * 4);

    for (int i = 0; i < msh->qud.len; ++i) {
        int b = i * 4;

        node[b + 0] = msh->qud.dat[i].vtx[1] + 1;
        node[b + 1] = msh->qud.dat[i].vtx[3] + 1;
        node[b + 2] = msh->qud.dat[i].vtx[2] + 1;
        node[b + 3] = msh->qud.dat[i].vtx[0] + 1;
    }

    cg_section_write(fi, bi, zi, "Elem", QUAD_4, 1, msh->qud.len, 0, node, &si);
    cg_close(fi);

    free(node);
    free(x);
    free(y);

    return 0;
}

static int umsh_exp_gns_c3d(struct umsh *msh, const char *fname)
{
    int fi;
    int bi;
    int zi;
    int ci;
    int si;

    int dim[3] = {msh->vtx.v3d.len, msh->hxd.len, 0};

    double *x = malloc(sizeof(double) * msh->vtx.v3d.len);
    double *y = malloc(sizeof(double) * msh->vtx.v3d.len);
    double *z = malloc(sizeof(double) * msh->vtx.v3d.len);

    for (int i = 0; i < msh->vtx.v3d.len; ++i) {
        x[i] = msh->vtx.v3d.dat[i].dat[0];
        y[i] = msh->vtx.v3d.dat[i].dat[1];
        z[i] = msh->vtx.v3d.dat[i].dat[2];
    }

    cg_open(fname, CG_MODE_WRITE, &fi);

    cg_base_write(fi, "Base", 3, 3, &bi);
    cg_zone_write(fi, bi, "Zone 1", dim, Unstructured, &zi);

    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateX", x, &ci);
    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateY", y, &ci);
    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateZ", z, &ci);

    int *node = malloc(sizeof(int *) * msh->hxd.len * 8);

    for (int i = 0; i < msh->hxd.len; ++i) {
        int b = i * 8;

        node[b + 0] = msh->hxd.dat[i].vtx[1] + 1;
        node[b + 1] = msh->hxd.dat[i].vtx[3] + 1;
        node[b + 2] = msh->hxd.dat[i].vtx[2] + 1;
        node[b + 3] = msh->hxd.dat[i].vtx[0] + 1;
        node[b + 4] = msh->hxd.dat[i].vtx[5] + 1;
        node[b + 5] = msh->hxd.dat[i].vtx[7] + 1;
        node[b + 6] = msh->hxd.dat[i].vtx[6] + 1;
        node[b + 7] = msh->hxd.dat[i].vtx[4] + 1;
    }

    cg_section_write(fi, bi, zi, "Elem", HEXA_8, 1, msh->hxd.len, 0, node, &si);
    cg_close(fi);

    free(node);
    free(x);
    free(y);
    free(z);

    return 0;
}

