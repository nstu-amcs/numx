#include <assert.h>

#include <numx/pde/msh.h>

cut_gen(vtx_cut, vtx, PUB);
cut_gen(seg_cut, seg, PUB);
cut_gen(qud_cut, qud, PUB);
cut_gen(hxd_cut, hxd, PUB);

int msh_imp(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    if (vtx_cut_new(&msh->vtx))
        return -1;

    if (seg_cut_new(&msh->seg))
        return -1;

    if (qud_cut_new(&msh->qud))
        return -1;

    if (hxd_cut_new(&msh->hxd))
        return -1;

    switch (msh->fmt) {
        case MSH_FMT_GRD:
            return msh_imp_grd(msh, dir, pfx);
        case MSH_FMT_GNS:
            return msh_imp_gns(msh, dir, pfx);
    }

    errno = ENOTSUP;
    return -1;
}

int msh_exp(struct msh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);

    switch (msh->fmt) {
        case MSH_FMT_GRD:
            return msh_exp_grd(msh, dir, pfx);
        case MSH_FMT_GNS:
            return msh_exp_gns(msh, dir, pfx);
    }

    errno = ENOTSUP;
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
