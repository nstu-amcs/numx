#include <assert.h>
#include <stdio.h>

#include <numx/pde/msh.h>

static int get_vtx(struct msh *msh, FILE *f);

static int get_qud_ems(struct msh *msh, FILE *f);
static int get_hxd_ems(struct msh *msh, FILE *f);

static int get_seg_bnd(struct msh *msh, FILE *f);
static int get_qud_bnd(struct msh *msh, FILE *f);

int msh_imp_grd(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    int r = 0;
    int n = 0;
    int e = 0;
    int b = 0;

    FILE *hdr = 0;
    FILE *vtx = 0;
    FILE *ems = 0;
    FILE *bnd = 0;

    char path[64];

    sprintf(path, "%s/%s.header", dir, pfx);

    if (!(hdr = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.nodes", dir, pfx);

    if (!(vtx = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.elements", dir, pfx);

    if (!(ems = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    sprintf(path, "%s/%s.boundary", dir, pfx);

    if (!(bnd = fopen(path, "r"))) {
        r = -1;
        goto end;
    }

    if (fscanf(hdr, "%d %d %d", &n, &e, &b) != 3) {
        r = -1;
        goto end;
    }

    switch (msh->sys) {
        case MSH_SYS_C2D:
            if ((r = vtx_cut_dev(&msh->vtx, n)))
                goto end;

            if ((r = qud_cut_dev(&msh->qud, e)))
                goto end;

            if ((r = seg_cut_dev(&msh->seg, b)))
                goto end;

            if ((r = get_seg_bnd(msh, bnd)))
                goto end;

            if ((r = get_qud_ems(msh, ems)))
                goto end;

            break;
        case MSH_SYS_C3D:
            if ((r = vtx_cut_dev(&msh->vtx, n)))
                goto end;

            if ((r = hxd_cut_dev(&msh->hxd, e)))
                goto end;

            if ((r = qud_cut_dev(&msh->qud, b)))
                goto end;

            if ((r = get_qud_bnd(msh, bnd)))
                goto end;

            if ((r = get_hxd_ems(msh, ems)))
                goto end;

            break;
    }

    if ((r = get_vtx(msh, vtx)))
        goto end;

    vtx_cut_shr(&msh->vtx);
    seg_cut_shr(&msh->seg);
    qud_cut_shr(&msh->qud);
    hxd_cut_shr(&msh->hxd);

end:
    fclose(hdr);
    fclose(vtx);
    fclose(ems);
    fclose(bnd);

    return r;
}

int msh_exp_grd(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    errno = ENOTSUP;

    return -1;
}

// clang-format off

static int get_vtx(struct msh *msh, FILE *f)
{
    struct vtx *vtx = msh->vtx.dat;

    for (int i = 0, j; i < msh->vtx.len; ++i)
        if (fscanf(f, "%d %d %lf %lf %lf", &j, &j, 
              &vtx[i].x, 
              &vtx[i].y, 
              &vtx[i].z) != 5)
            return -1;

    return 0;
}

static int get_qud_ems(struct msh *msh, FILE *f)
{
    struct qud *qud = msh->qud.dat;

    for (int i = 0, j; i < msh->qud.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d", &j, &qud[i].pid, &j, 
              &qud[i].vtx[0], 
              &qud[i].vtx[1], 
              &qud[i].vtx[2],
              &qud[i].vtx[3]) != 7)
            return -1;

        qud[i].pid -= 1;
        qud[i].vtx[0] -= 1;
        qud[i].vtx[1] -= 1;
        qud[i].vtx[2] -= 1;
        qud[i].vtx[3] -= 1;
    }

    return 0;
}

static int get_hxd_ems(struct msh *msh, FILE *f)
{
    struct hxd *hxd = msh->hxd.dat;

    for (int i = 0, j; i < msh->hxd.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d", &j, &hxd[i].pid, &j, 
              &hxd[i].vtx[0], 
              &hxd[i].vtx[1],
              &hxd[i].vtx[3],
              &hxd[i].vtx[2], 
              &hxd[i].vtx[4], 
              &hxd[i].vtx[5], 
              &hxd[i].vtx[7], 
              &hxd[i].vtx[6]) != 11)
            return -1;

        hxd[i].pid -= 1;
        hxd[i].vtx[0] -= 1;
        hxd[i].vtx[1] -= 1;
        hxd[i].vtx[2] -= 1;
        hxd[i].vtx[3] -= 1;
        hxd[i].vtx[4] -= 1;
        hxd[i].vtx[5] -= 1;
        hxd[i].vtx[6] -= 1;
        hxd[i].vtx[7] -= 1;
    }

    return 0;
}

static int get_seg_bnd(struct msh *msh, FILE *f)
{
    struct seg *seg = msh->seg.dat;

    for (int i = 0, j; i < msh->seg.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d", &j, &seg[i].pid, &j, &j, &j, 
              &seg[i].vtx[0], 
              &seg[i].vtx[1]) != 7)
            return -1;

        seg[i].pid -= 1;
        seg[i].vtx[0] -= 1;
        seg[i].vtx[1] -= 1;
    }

    return 0;
}

static int get_qud_bnd(struct msh *msh, FILE *f)
{
    struct qud *qud = msh->qud.dat;

    for (int i = 0, j; i < msh->qud.len; ++i) {
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &j, &j, 
              &qud[i].vtx[0], 
              &qud[i].vtx[1],
              &qud[i].vtx[2], 
              &qud[i].vtx[3]) != 9)
            return -1;

        qud[i].pid -= 1;
        qud[i].vtx[0] -= 1;
        qud[i].vtx[1] -= 1;
        qud[i].vtx[2] -= 1;
        qud[i].vtx[3] -= 1;
    }

    return 0;
}
