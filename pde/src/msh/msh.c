#include <assert.h>

#include <numx/pde/msh.h>

cut_gen(seg_cut, seg, PUB);
cut_gen(qud_cut, qud, PUB);
cut_gen(hxd_cut, hxd, PUB);

int msh_new(struct msh *msh)
{
    assert(msh);

    if (vec_cut_new(&msh->vtx))
        goto err;

    if (seg_cut_new(&msh->seg))
        goto err;

    if (qud_cut_new(&msh->qud))
        goto err;

    if (hxd_cut_new(&msh->hxd))
        goto err;

    return 0;

err:
    msh_cls(msh);
    return -1;
}

int msh_cls(struct msh *msh)
{
    assert(msh);

    for (int i = 0; i < msh->vtx.len; ++i)
        vec_cls(&msh->vtx.dat[i]);

    vec_cut_cls(&msh->vtx);
    seg_cut_cls(&msh->seg);
    qud_cut_cls(&msh->qud);
    hxd_cut_cls(&msh->hxd);

    return 0;
}

int msh_qud_nrm(struct msh *msh, struct qud *qud, struct vec *nrm)
{
    struct vec *vtx = msh->vtx.dat;

    struct vec *a = &vtx[qud->vtx[0]];
    struct vec *b = &vtx[qud->vtx[1]];
    struct vec *c = &vtx[qud->vtx[2]];

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
