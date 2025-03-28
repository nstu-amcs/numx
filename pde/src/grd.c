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

int msh_new(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);

    int r = 0;

    if ((r = vtx_cut_new(&msh->vtx)))
        return r;

    if ((r = seg_cut_new(&msh->seg)))
        return r;

    if ((r = qud_cut_new(&msh->qud)))
        return r;

    if ((r = hxd_cut_new(&msh->hxd)))
        return r;

    FILE *hdr = 0;
    FILE *vtx = 0;
    FILE *ems = 0;
    FILE *bnd = 0;

    int  nc, ec, bc;
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

    if (fscanf(hdr, "%d %d %d", &nc, &ec, &bc) != 3) {
        r = -1;
        goto end;
    }

    switch (msh->sys) {
        case MSH_C2D:
            if ((r = vtx_cut_dev(&msh->vtx, nc)))
                goto end;

            if ((r = qud_cut_dev(&msh->qud, ec)))
                goto end;

            if ((r = seg_cut_dev(&msh->seg, bc)))
                goto end;

            break;
        case MSH_C3D:
            if ((r = vtx_cut_dev(&msh->vtx, nc)))
                goto end;

            if ((r = hxd_cut_dev(&msh->hxd, ec)))
                goto end;

            if ((r = qud_cut_dev(&msh->qud, bc)))
                goto end;

            break;
    }

    if ((r = get_vtx(msh, vtx)))
        goto end;

    if ((r = get_ems(msh, ems)))
        goto end;

    if ((r = get_bnd(msh, bnd)))
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

    if (r) {
        vtx_cut_cls(&msh->vtx);
        seg_cut_cls(&msh->seg);
        qud_cut_cls(&msh->qud);
        hxd_cut_cls(&msh->hxd);
    }

    return r;
}

static int get_vtx(struct msh *msh, FILE *f)
{
    struct vtx *vtx = msh->vtx.dat;

    for (int i = 0, j; i < msh->vtx.len; ++i)
        if (fscanf(f, "%d %d %lf %lf %lf", &j, &j, &vtx[i].x, &vtx[i].y, &vtx[i].z) != 5)
            return -1;

    return 0;
}

static int get_ems(struct msh *msh, FILE *f)
{
    struct qud *qud = msh->qud.dat;
    struct hxd *hxd = msh->hxd.dat;

    switch (msh->sys) {
        case MSH_C2D:
            for (int i = 0, j; i < msh->qud.len; ++i) {
                if (fscanf(f, "%d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &qud[i].vtx[0], &qud[i].vtx[1],
                        &qud[i].vtx[2], &qud[i].vtx[3]) != 7)
                    return -1;

                qud[i].pid -= 1;
                qud[i].vtx[0] -= 1;
                qud[i].vtx[1] -= 1;
                qud[i].vtx[2] -= 1;
                qud[i].vtx[3] -= 1;
            }

            break;
        case MSH_C3D:
            for (int i = 0, j; i < msh->hxd.len; ++i) {
                if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d", &j, &hxd[i].pid, &j, &hxd[i].vtx[0], &hxd[i].vtx[1],
                        &hxd[i].vtx[2], &hxd[i].vtx[3], &hxd[i].vtx[4], &hxd[i].vtx[5], &hxd[i].vtx[6],
                        &hxd[i].vtx[7]) != 11)
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

            break;
    }

    return 0;
}

static int get_bnd(struct msh *msh, FILE *f)
{
    struct seg *seg = msh->seg.dat;
    struct qud *qud = msh->qud.dat;

    switch (msh->sys) {
        case MSH_C2D:
            for (int i = 0, j; i < msh->seg.len; ++i) {
                if (fscanf(f, "%d %d %d %d %d %d %d", &j, &seg[i].pid, &j, &j, &j, &seg[i].vtx[0], &seg[i].vtx[1]) != 7)
                    return -1;

                seg[i].pid -= 1;
                seg[i].vtx[0] -= 1;
                seg[i].vtx[1] -= 1;
            }

            break;
        case MSH_C3D:
            for (int i = 0, j; i < msh->qud.len; ++i) {
                if (fscanf(f, "%d %d %d %d %d %d %d %d %d", &j, &qud[i].pid, &j, &j, &j, &qud[i].vtx[0], &qud[i].vtx[1],
                        &qud[i].vtx[2], &qud[i].vtx[3]) != 9)
                    return -1;

                qud[i].pid -= 1;
                qud[i].vtx[0] -= 1;
                qud[i].vtx[1] -= 1;
                qud[i].vtx[2] -= 1;
                qud[i].vtx[3] -= 1;
            }

            break;
    }

    return 0;
}
