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
        return -1;

    if (seg_cut_new(&msh->seg))
        return -1;

    if (qud_cut_new(&msh->qud))
        return -1;

    if (hxd_cut_new(&msh->hxd))
        return -1;

    return 0;
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
