#include <assert.h>
#include <dlfcn.h>

#include <numx/pde/cnd.h>
#include <numx/pde/sim.h>

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
    sim->ops.usr = 0;
    sim->ops.exp.mod = SIM_EXP_GNS;
    sim->ops.tdd.num = 0;
    sim->ops.tdd.hop = 0;

    if (!(sim->msh = malloc(sizeof(struct msh))))
        goto err;

    if (msh_new(sim->msh))
        goto err;

    if (mat_cut_new(&sim->mat))
        goto err;

    if (val_cut_new(&sim->src))
        goto err;

    if (obj_cut_new(&sim->obj))
        goto err;

    if (bnd_cut_new(&sim->bnd))
        goto err;

    if (cnd_ini_cut_new(&sim->cnd_ini))
        goto err;

    if (cnd_bnd_cut_new(&sim->cnd_bnd))
        goto err;

    return 0;

err:
    sim_cls(sim);

    return -1;
}

int sim_cls(struct sim *sim)
{
    assert(sim);

    msh_cls(sim->msh);

    free(sim->msh);

    mat_cut_cls(&sim->mat);
    val_cut_cls(&sim->src);
    obj_cut_cls(&sim->obj);
    bnd_cut_cls(&sim->bnd);

    cnd_ini_cut_cls(&sim->cnd_ini);
    cnd_bnd_cut_cls(&sim->cnd_bnd);

    dlclose(sim->ops.usr);

    return 0;
}

int sim_run(struct sim *sim)
{
    assert(sim);
    return sim->slv->exe(sim);
}
