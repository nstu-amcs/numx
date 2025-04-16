#include <assert.h>

#include "fem_lin.h"

static int fem_lin_std_slv(struct sim *sim);
static int fem_lin_non_slv(struct sim *sim);
static int fem_lin_hmc_slv(struct sim *sim);

int fem_lin_ell_slv(struct sim *sim)
{
    assert(sim);

    int r = 0;

    struct fem *fem = (struct fem *)sim->slv;
    struct vec  wgt;

    if ((r = vec_new(&wgt, fem->prv.vec.n)))
        goto end;

    sim->rt.ti = 0;
    sim->rt.tv = 0;
    sim->rt.wgt = &wgt;

    switch (fem->ops.mod) {
        case FEM_STD:
            if ((r = fem_lin_std_slv(sim)))
                goto end;

            break;
        case FEM_NON:
            if ((r = fem_lin_non_slv(sim)))
                goto end;

            break;
        case FEM_HMC:
            if ((r = fem_lin_hmc_slv(sim)))
                goto end;

            break;
    }

    if (fem->slv.itr.run)
        fem->slv.itr.run(fem->slv.itr.ctx, sim);

    if (sim->ops.exp.put.run)
        sim->ops.exp.put.run(sim->ops.exp.put.ctx, sim);

end:
    vec_cls(&wgt);
    return r;
}

static int fem_lin_std_slv(struct sim *sim)
{
    struct fem *fem = (struct fem *)sim->slv;

    if (fem_lin_ell_asm(sim))
        return -1;

    switch (fem->ops.iss.mod) {
        case ISS_BCG:
            if (iss_bcg_slv(&fem->prv.ell, sim->rt.wgt, &fem->prv.vec, &fem->ops.iss.ops.bcg))
                return -1;

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int fem_lin_non_slv(struct sim *sim)
{
    int r = 0;

    struct fem  *fem = (struct fem *)sim->slv;
    struct smtx *mtx = &fem->prv.ell;
    struct vec  *vec = &fem->prv.vec;
    struct vec   tmp;

    if ((r = vec_new(&tmp, vec->n)))
        goto end;

    int    max = fem->ops.non.max;
    double err = fem->ops.non.err;
    double cur = 1;
    double tn = 0;
    double bn = 0;

    for (int i = 0; i < max && cur > err; ++i) {
        if (i > 0)
            switch (fem->ops.iss.mod) {
                case ISS_BCG:
                    if ((r = iss_bcg_slv(mtx, sim->rt.wgt, vec, &fem->ops.iss.ops.bcg)))
                        goto end;

                    break;
                default:
                    errno = ENOTSUP;
                    r = -1;
                    goto end;
            }

        if (fem_lin_ell_asm(sim))
            return -1;

        mtx_vmlt(mtx, sim->rt.wgt, &tmp);
        vec_cmb(&tmp, vec, &tmp, -1);
        vec_nrm(&tmp, &tn);
        vec_nrm(vec, &bn);

        cur = tn / bn;

        fem->ops.non.rt.itr = i;
        fem->ops.non.rt.err = cur;

        if (fem->ops.non.itr.run)
            fem->ops.non.itr.run(fem->ops.non.itr.ctx, sim);
    }

end:
    vec_cls(&tmp);

    return r;
}

static int fem_lin_hmc_slv(struct sim *sim)
{
    assert(sim);

    return 0;
}
