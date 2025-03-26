#include <assert.h>
#include <stdio.h>

#include <numx/com/dev/cut.h>
#include <numx/pde/msh.h>

cut_gen(vtx_cut, vtx, PUB);
cut_gen(seg_cut, seg, PUB);
cut_gen(qud_cut, qud, PUB);
cut_gen(hxd_cut, hxd, PUB);

static int get_vtx(struct msh *msh, FILE *f);
static int get_ems(struct msh *msh, FILE *f);
static int get_bnd(struct msh *msh, FILE *f);

int msh_new(struct msh *msh, const char *dir)
{
    assert(msh);
    assert(dir);

    int r = 0;

    vtx_cut_new(&msh->vtx);
    seg_cut_new(&msh->seg);
    qud_cut_new(&msh->qud);
    hxd_cut_new(&msh->hxd);

    FILE *hdr = 0;
    FILE *vtx = 0;
    FILE *ems = 0;
    FILE *bnd = 0;

    int  nc, ec, bc;
    char path[64];

    sprintf(path, "%s/mesh.header", dir);

    if (!(hdr = fopen(path, "r")))
        goto end;

    sprintf(path, "%s/mesh.nodes", dir);

    if (!(vtx = fopen(path, "r")))
        goto end;

    sprintf(path, "%s/mesh.elements", dir);

    if (!(ems = fopen(path, "r")))
        goto end;

    sprintf(path, "%s/mesh.boundary", dir);

    if (!(bnd = fopen(path, "r")))
        goto end;

    if (fscanf(hdr, "%d %d %d", &nc, &ec, &bc) != 3)
        goto end;

    switch (msh->sys) {
        case MSH_C2D:
            if (vtx_cut_dev(&msh->vtx, nc))
                goto end;

            if (qud_cut_dev(&msh->qud, ec))
                goto end;

            if (seg_cut_dev(&msh->seg, bc))
                goto end;

            break;
        case MSH_C3D:
            if (vtx_cut_dev(&msh->vtx, nc))
                goto end;

            if (hxd_cut_dev(&msh->hxd, ec))
                goto end;

            if (qud_cut_dev(&msh->qud, bc))
                goto end;

            break;
    }

    if ((r = get_vtx(msh, vtx)))
        goto end;

    if ((r = get_ems(msh, ems)))
        goto end;

    if ((r = get_bnd(msh, bnd)))
        goto end;

end:
    fclose(hdr);
    fclose(vtx);
    fclose(ems);
    fclose(bnd);

    return r;
}

static int get_vtx(struct msh *msh, FILE *f)
{
    vtx *vtx = msh->vtx.dat;

    for (int i = 0, j; i < msh->vtx.len; ++i)
        if (fscanf(f, "%d %d %lf %lf %lf", &j, &j, &vtx[i].x, &vtx[i].y, &vtx[i].z) != 5)
            return -1;

    return 0;
}

static int get_ems(struct msh *msh, FILE *f)
{
    qud *qud = msh->qud.dat;
    hxd *hxd = msh->hxd.dat;

    switch (msh->sys) {
        case MSH_C2D:
            for (int i = 0, j; i < msh->qud.len; ++i)
                if (fscanf(f, "%d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &qud[i].vtx[0], &qud[i].vtx[1],
                        &qud[i].vtx[2], &qud[i].vtx[3]) != 7)
                    return -1;

            break;
        case MSH_C3D:
            for (int i = 0, j; i < msh->hxd.len; ++i)
                if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d", &j, &hxd[i].pid, &j, &hxd[i].vtx[0], &hxd[i].vtx[1],
                        &hxd[i].vtx[2], &hxd[i].vtx[3], &hxd[i].vtx[4], &hxd[i].vtx[5], &hxd[i].vtx[6],
                        &hxd[i].vtx[7]) != 11)
                    return -1;

            break;
    }

    return 0;
}

static int get_bnd(struct msh *msh, FILE *f)
{
    seg *seg = msh->seg.dat;
    qud *qud = msh->qud.dat;

    switch (msh->sys) {
        case MSH_C2D:
            for (int i = 0, j; i < msh->seg.len; ++i)
                if (fscanf(f, "%d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &j, &j, &seg[i].vtx[0], &seg[i].vtx[1]) != 7)
                    return -1;

            break;
        case MSH_C3D:
            for (int i = 0, j; i < msh->hxd.len; ++i)
                if (fscanf(f, "%d %d %d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &j, &j, &qud[i].vtx[0], &qud[i].vtx[1],
                        &qud[i].vtx[2], &qud[i].vtx[3]) != 9)
                    return -1;

            break;
    }

    return 0;
}
