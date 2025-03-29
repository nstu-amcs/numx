#include <assert.h>

#include <numx/pde/sim.h>
#include <numx/pde/val.h>

cut_gen(val_cut, val, PUB);
cut_gen(mat_cut, mat, PUB);
cut_gen(obj_cut, obj, PUB);
cut_gen(bnd_cut, bnd, PUB);

cut_gen(cnd_bnd_cut, cnd_bnd, PUB);
cut_gen(cnd_ini_cut, cnd_ini, PUB);

int sim_new(struct sim *sim)
{
    assert(sim);

    sim->mod = SIM_ELL;

    sim->ops.ell.ops.usr = 0;
    sim->ops.ell.ops.exp.mod = SIM_EXP_GNS;

    sim->ops.ell.ops.non.max = 100;
    sim->ops.ell.ops.non.err = 1e10;
    sim->ops.ell.ops.non.rlx = 1;

    sim->ops.ell.ops.iss.mod = ISS_BCG;
    sim->ops.ell.ops.iss.ops.bcg.con.sm = 0;
    sim->ops.ell.ops.iss.ops.bcg.ops.max = 500;
    sim->ops.ell.ops.iss.ops.bcg.ops.err = 1e10;
    sim->ops.ell.ops.iss.ops.bcg.ops.itr.run = 0;

    sim->ops.pbc.num = 1;
    sim->ops.pbc.hop = 0;

    if (!(sim->msh = malloc(sizeof(struct msh))))
        return -1;

    if (mat_cut_new(&sim->mat))
        return -1;

    if (val_cut_new(&sim->src))
        return -1;

    if (obj_cut_new(&sim->obj))
        return -1;

    if (bnd_cut_new(&sim->bnd))
        return -1;

    if (cnd_ini_cut_new(&sim->cnd_ini))
        return -1;

    if (cnd_bnd_cut_new(&sim->cnd_bnd))
        return -1;

    return 0;
}

int sim_cls(struct sim *sim)
{
    assert(sim);

    msh_cls(sim->msh);

    mat_cut_cls(&sim->mat);
    val_cut_cls(&sim->src);
    obj_cut_cls(&sim->obj);
    bnd_cut_cls(&sim->bnd);

    cnd_ini_cut_cls(&sim->cnd_ini);
    cnd_bnd_cut_cls(&sim->cnd_bnd);

    free(sim->msh);

    return 0;
}
