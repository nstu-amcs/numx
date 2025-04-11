#include <assert.h>
#include <cgnslib.h>
#include <stdio.h>

#include <numx/pde/sim.h>

int sim_exp_ini_gns(void *ctx, struct sim *sim)
{
    (void)ctx;

    assert(sim);
    return msh_exp_gns(sim->msh, sim->ops.exp.dir, sim->ops.exp.pfx);
}

int sim_exp_put_gns(void *ctx, struct sim *sim)
{
    (void)ctx;

    assert(sim);

    char fname[256];
    char sname[64];

    int fi;
    int si;
    int ii;

    sprintf(fname, "%s/%s.cgns", sim->ops.exp.dir, sim->ops.exp.pfx);
    sprintf(sname, "FlowSolution %d", sim->rt.tm);

    cg_open(fname, CG_MODE_MODIFY, &fi);
    cg_sol_write(fi, 1, 1, sname, Vertex, &si);
    cg_field_write(fi, 1, 1, si, RealDouble, "Temperature", sim->rt.wgt->dat, &ii);
    cg_close(fi);

    return 0;
}
