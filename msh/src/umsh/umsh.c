#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>

cut_gen(v2d_cut, v2d, PUB);
cut_gen(v3d_cut, v3d, PUB);
cut_gen(seg_cut, seg, PUB);
cut_gen(qud_cut, qud, PUB);
cut_gen(hxd_cut, hxd, PUB);

int umsh_new(struct umsh *msh)
{
    assert(msh);

    msh->type = MSH_C2D;

    if (v2d_cut_new(&msh->vtx.v2d))
        goto err;

    if (seg_cut_new(&msh->seg))
        goto err;

    if (qud_cut_new(&msh->qud))
        goto err;

    if (hxd_cut_new(&msh->hxd))
        goto err;

    return 0;

err:
    umsh_cls(msh);
    return -1;
}

int umsh_cls(struct umsh *msh)
{
    assert(msh);

    switch (msh->type) {
        case MSH_C2D:
            v2d_cut_cls(&msh->vtx.v2d);
            break;
        case MSH_C3D:
            v3d_cut_cls(&msh->vtx.v3d);
            break;
    }

    seg_cut_cls(&msh->seg);
    qud_cut_cls(&msh->qud);
    hxd_cut_cls(&msh->hxd);

    return 0;
}

int umsh_seg_nrm(struct umsh *msh, struct seg *seg, struct vec *nrm)
{
    struct v2d *vtx = msh->vtx.v2d.dat;

    struct v2d *a = &vtx[seg->vtx[0]];
    struct v2d *b = &vtx[seg->vtx[1]];

    memset(nrm->dat, 0, sizeof(double) * 2);

    if (a->dat[0] == b->dat[0]) {
        nrm->dat[0] = 1;
    } else {
        nrm->dat[1] = 1;
    }

    return 0;
}

int umsh_qud_nrm(struct umsh *msh, struct qud *qud, struct vec *nrm)
{
    struct v3d *vtx = msh->vtx.v3d.dat;

    struct v3d *a = &vtx[qud->vtx[0]];
    struct v3d *b = &vtx[qud->vtx[1]];
    struct v3d *c = &vtx[qud->vtx[2]];

    memset(nrm->dat, 0, sizeof(double) * 3);

    if (a->dat[0] == c->dat[0] && b->dat[0] == c->dat[0]) {
        nrm->dat[0] = 1;
        return 0;
    }

    if (a->dat[1] == c->dat[1] && b->dat[1] == c->dat[1]) {
        nrm->dat[1] = 1;
        return 0;
    }

    nrm->dat[2] = 1;
    return 0;
}

int umsh_qud_vtx_loc(struct qud *qud, int gv)
{
    for (int i = 0; i < 4; ++i) {
        if (qud->vtx[i] == gv) {
            return i;
        }
    }

    return -1;
}

static int qud_bnd_cmb[4][2] = {
    {0, 1},
    {2, 3},
    {0, 2},
    {1, 3},
};

int umsh_seg_srh(struct umsh *msh, seg_srh_fun fun, int pid)
{
    vtx_ptr vtx[4];

    for (int i = 0; i < msh->qud.len; ++i) {
        struct qud *q = &msh->qud.dat[i];

        switch (msh->type) {
            case MSH_C2D:
                vtx[0].v2d = &msh->vtx.v2d.dat[q->vtx[0]];
                vtx[1].v2d = &msh->vtx.v2d.dat[q->vtx[1]];
                vtx[2].v2d = &msh->vtx.v2d.dat[q->vtx[2]];
                vtx[3].v2d = &msh->vtx.v2d.dat[q->vtx[3]];
                break;
            case MSH_C3D:
                vtx[0].v3d = &msh->vtx.v3d.dat[q->vtx[0]];
                vtx[1].v3d = &msh->vtx.v3d.dat[q->vtx[1]];
                vtx[2].v3d = &msh->vtx.v3d.dat[q->vtx[2]];
                vtx[3].v3d = &msh->vtx.v3d.dat[q->vtx[3]];
                break;
        }

        struct seg s = {
            .pid = pid,
            .qud = i,
        };

        for (int j = 0; j < 4; ++j) {
            int a = qud_bnd_cmb[j][0];
            int b = qud_bnd_cmb[j][1];

            s.vtx[0] = q->vtx[a];
            s.vtx[1] = q->vtx[b];

            if (fun(msh, &s)) {
                seg_cut_add(&msh->seg, s);
            }
        }
    }

    return 0;
}

static int umsh_c2d_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx);
static int umsh_c3d_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx);

int umsh_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx)
{
    switch (msh->type) {
        case MSH_C2D:
            return umsh_c2d_vtx_qud_lup(msh, vtx);
        case MSH_C3D:
            return umsh_c3d_vtx_qud_lup(msh, vtx);
    }

    return -1;
}

static int umsh_c2d_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx)
{
    static const double tol = 1e-10;

    for (int qi = 0; qi < msh->qud.len; ++qi) {
        struct qud *qud = &msh->qud.dat[qi];
        struct v2d *a = &msh->vtx.v2d.dat[qud->vtx[0]];
        struct v2d *b = &msh->vtx.v2d.dat[qud->vtx[3]];
        struct v2d *v = vtx.v2d;

        double x0 = a->dat[0];
        double y0 = a->dat[1];
        double x1 = b->dat[0];
        double y1 = b->dat[1];

        double x = v->dat[0];
        double y = v->dat[1];

        if (less(x, x0, tol) || less(x1, x, tol)) {
            continue;
        }

        if (less(y, y0, tol) || less(y1, y, tol)) {
            continue;
        }

        return qi;
    }

    return -1;
}

static int umsh_c3d_vtx_qud_lup(struct umsh *msh, union vtx_ptr vtx)
{
    return -1;
}