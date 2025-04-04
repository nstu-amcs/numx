#include <assert.h>
#include <cgnslib.h>
#include <stdarg.h>
#include <stdio.h>

#include <numx/pde/sim.h>

int sim_exp_ini_gns(void *ctx, int n, ...)
{
    assert(ctx);
    assert(n == 0);

    struct sim *sim = (struct sim *)ctx;

    return msh_exp_gns(sim->msh, sim->ops.exp.dir, sim->ops.exp.pfx);
}

int sim_exp_put_gns(void *ctx, int n, ...)
{
    assert(ctx);
    assert(n == 2);

    va_list arg;
    va_start(arg, n);

    int i = va_arg(arg, int);

    struct vec *wgt = va_arg(arg, struct vec *);
    struct sim *sim = (struct sim *)ctx;

    va_end(arg);

    char fname[256];
    char sname[64];

    int fi;
    int si;
    int ii;

    sprintf(fname, "%s/%s.cgns", sim->ops.exp.dir, sim->ops.exp.pfx);
    sprintf(sname, "FlowSolution %d", i);

    cg_open(fname, CG_MODE_MODIFY, &fi);
    cg_sol_write(fi, 1, 1, sname, Vertex, &si);
    cg_field_write(fi, 1, 1, si, RealDouble, "Temperature", wgt->dat, &ii);
    cg_close(fi);

    return 0;
}
