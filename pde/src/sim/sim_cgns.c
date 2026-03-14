#include </usr/include/cgns_io.h>
#include </usr/include/cgnslib.h>
#include <assert.h>
#include <stdio.h>

#include <numx/pde/sim.h>

int sim_exp_cgns_ini(struct sim *sim)
{
    if (sim->mod == SIM_HMC)
        return 0;

    if (umsh_exp_cgns(sim->msh, sim->ops.exp.dir, sim->ops.exp.pfx))
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

        sprintf(ss, "FlowSolution-%d-%d", 0, i);
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

int sim_exp_cgns_put_v(struct sim *sim, const char *name, int i, struct vec *s)
{
    if (sim->mod == SIM_HMC)
        return 0;

    char gname[256];
    char fname[256];
    char sname[64];

    int fi;
    int bi;
    int zi;
    int si;
    int ii;
    int io;

    double root_id;
    double zone_id;
    double link_id;

    int dim[3] = {sim->msh->vtx.v2d.len, sim->msh->qud.len, 0};

    sprintf(gname, "%s/%s.cgns", sim->ops.exp.dir, sim->ops.exp.pfx);
    sprintf(fname, "%s/%s-%d.cgns", sim->ops.exp.dir, sim->ops.exp.pfx, i);
    sprintf(sname, "FlowSolution-%d", sim->slv->run.ti);

    if (cg_open(fname, CG_MODE_WRITE, &fi) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    if (cg_base_write(fi, "Base", 2, 2, &bi) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    if (cg_zone_write(fi, bi, "Zone1", dim, Unstructured, &zi) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    cg_get_cgio(fi, &io);
    cgio_get_root_id(io, &root_id);
    cgio_get_node_id(io, root_id, "/Base/Zone1", &zone_id);

    if (cgio_create_link(
            io, zone_id, "GridCoordinates", gname, "/Base/Zone1/GridCoordinates", &link_id)) {
        cg_error_print();
        return -EIO;
    }

    if (cgio_create_link(
            io, zone_id, "Elem", gname, "/Base/Zone1/Elem", &link_id)) {
        cg_error_print();
        return -EIO;
    }

    if (cg_sol_write(fi, 1, 1, sname, Vertex, &si) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    if (cg_field_write(fi, 1, 1, si, RealDouble, name, s->dat, &ii) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    if (cg_close(fi) != CG_OK) {
        cg_error_print();
        return -EIO;
    }

    return 0;
}

int sim_exp_cgns_put_c(struct sim *sim, const char *name, int i, struct vec *s)
{
    if (sim->mod == SIM_HMC)
        return 0;

    char fname[256];
    char sname[64];

    int fi;
    int si;
    int ii;

    sprintf(fname, "%s/%s-%d.cgns", sim->ops.exp.dir, sim->ops.exp.pfx, i);
    sprintf(sname, "FlowSolutionCell-%d", sim->slv->run.ti);

    cg_open(fname, CG_MODE_MODIFY, &fi);
    cg_sol_write(fi, 1, 1, sname, CellCenter, &si);
    cg_field_write(fi, 1, 1, si, RealDouble, name, s->dat, &ii);
    cg_close(fi);

    return 0;
}
