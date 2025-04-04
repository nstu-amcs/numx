#include <assert.h>

#include <numx/pde/msh.h>

cut_gen(vtx_cut, vtx, PUB);
cut_gen(seg_cut, seg, PUB);
cut_gen(qud_cut, qud, PUB);
cut_gen(hxd_cut, hxd, PUB);

int msh_new(struct msh *msh)
{
    assert(msh);

    if (vtx_cut_new(&msh->vtx))
        goto new_err;

    if (seg_cut_new(&msh->seg))
        goto new_err;

    if (qud_cut_new(&msh->qud))
        goto new_err;

    if (hxd_cut_new(&msh->hxd))
        goto new_err;

    return 0;

new_err:
    msh_cls(msh);

    return -1;
}

int msh_cls(struct msh *msh)
{
    assert(msh);

    vtx_cut_cls(&msh->vtx);
    seg_cut_cls(&msh->seg);
    qud_cut_cls(&msh->qud);
    hxd_cut_cls(&msh->hxd);

    return 0;
}

int msh_qud_nrm(struct msh *msh, int q)
{
    struct qud *qud = &msh->qud.dat[q];
    struct vtx *vtx = msh->vtx.dat;

    struct vtx *a = &vtx[qud->vtx[0]];
    struct vtx *b = &vtx[qud->vtx[1]];
    struct vtx *c = &vtx[qud->vtx[2]];

    if (a->x == c->x && b->x == c->x)
        return 0;

    if (a->y == c->y && b->y == c->y)
        return 1;

    return 2;
}
