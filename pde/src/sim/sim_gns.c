#include <assert.h>
#include <cgnslib.h>
#include <stdio.h>

#include <numx/pde/sim.h>

int sim_exp_gns_ini(struct sim *sim)
{
    if (sim->mod == SIM_HMC)
        return 0;

    if (msh_exp_gns(sim->msh, sim->ops.exp.dir, sim->ops.exp.pfx))
        return -1;

    if (sim->eqn == SIM_ELL)
        return 0;

    int num = sim->ops.tdd.num;
    int hop = sim->ops.tdd.hop;

    char fname[256];

    int fi;
    int nu = num + 1;

    double *tm = malloc(sizeof(double) * (num + 1));
    char   *ss = malloc(sizeof(char) * 32);
    char   *sn = malloc(sizeof(char) * ((num + 1) * 32 + 1));

    for (int i = 0; i <= num; ++i) {
        tm[i] = hop * i;

        sprintf(ss, "FlowSolution-%d", i);
        sprintf(sn + (32 * i), "%-32s", ss);
    }

    int info[2] = {32, nu};

    sprintf(fname, "%s/%s.cgns", sim->ops.exp.dir, sim->ops.exp.pfx);

    cg_open(fname, CG_MODE_MODIFY, &fi);

    cg_biter_write(fi, 1, "TimeIterValues", num + 1);
    cg_goto(fi, 1, "BaseIterativeData_t", 1, "end");
    cg_array_write("TimeValues", RealDouble, 1, &nu, tm);

    cg_ziter_write(fi, 1, 1, "ZoneIterativeData");
    cg_goto(fi, 1, "Zone_t", 1, "ZoneIterativeData_t", 1, "end");
    cg_array_write("FlowSolutionPointers", Character, 2, info, sn);

    cg_simulation_type_write(fi, 1, TimeAccurate);

    cg_close(fi);

    free(tm);
    free(ss);
    free(sn);

    return 0;
}

int sim_exp_gns_put(struct sim *sim)
{
    if (sim->mod == SIM_HMC)
        return 0;

    char fname[256];
    char sname[64];

    int fi;
    int si;
    int ii;

    sprintf(fname, "%s/%s.cgns", sim->ops.exp.dir, sim->ops.exp.pfx);
    sprintf(sname, "FlowSolution-%d", sim->slv->run.ti);

    cg_open(fname, CG_MODE_MODIFY, &fi);
    cg_sol_write(fi, 1, 1, sname, Vertex, &si);
    cg_field_write(fi, 1, 1, si, RealDouble, "Temperature",
        sim->slv->run.wgt[0]->dat, &ii);
    cg_close(fi);

    return 0;
}
