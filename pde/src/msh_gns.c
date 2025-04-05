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
