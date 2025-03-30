#include <assert.h>
#include <cgnslib.h>
#include <stdio.h>
#include <stdlib.h>

#include <numx/pde/msh.h>

int msh_exp_gns(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    char fname[256];

    int fi;
    int bi;
    int zi;
    int ci;
    int si;

    int dim[3] = {msh->vtx.len, msh->hxd.len, 0};

    double *x = malloc(sizeof(double) * msh->vtx.len);
    double *y = malloc(sizeof(double) * msh->vtx.len);
    double *z = malloc(sizeof(double) * msh->vtx.len);

    for (int i = 0; i < msh->vtx.len; ++i) {
        x[i] = msh->vtx.dat[i].x;
        y[i] = msh->vtx.dat[i].y;
        z[i] = msh->vtx.dat[i].z;
    }

    sprintf(fname, "%s/%s.cgns", dir, pfx);
    cg_open(fname, CG_MODE_WRITE, &fi);

    cg_base_write(fi, "Base", 3, 3, &bi);
    cg_zone_write(fi, bi, "Zone 1", dim, Unstructured, &zi);

    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateX", x, &ci);
    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateY", y, &ci);
    cg_coord_write(fi, bi, zi, RealDouble, "CoordinateZ", z, &ci);

    char name[32];
    int  node[8];

    for (int i = 0; i < msh->hxd.len; ++i) {
        sprintf(name, "Elem %d", i);

        node[0] = msh->hxd.dat[i].vtx[1];
        node[1] = msh->hxd.dat[i].vtx[3];
        node[2] = msh->hxd.dat[i].vtx[2];
        node[3] = msh->hxd.dat[i].vtx[0];
        node[4] = msh->hxd.dat[i].vtx[5];
        node[5] = msh->hxd.dat[i].vtx[7];
        node[6] = msh->hxd.dat[i].vtx[6];
        node[7] = msh->hxd.dat[i].vtx[4];

        cg_section_write(fi, bi, zi, name, HEXA_8, 1, 1, 0, node, &si);
    }

    cg_close(fi);

    return 0;
}
