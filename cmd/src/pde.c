#include <dlfcn.h>
#include <math.h>
#include <stdio.h>

#include <numx/cmd.h>
#include <numx/pde/fem.h>
#include <numx/pde/sim.h>

#define XS 0.1
#define YS 0.1
#define ZS 0.01

static void man()
{
    printf("usage: numx pde {filename.sif}\n");
}

struct dat
{
    double (*tgt)(struct sim_fun_ctx *, struct vec *);

    FILE *fs;
    FILE *fn;

    struct
    {
        int    itr;
        double err;
    } non;
};

void on_slv(void *ctx, struct sim *sim)
{
    struct dat *dat = (struct dat *)ctx;
    struct vec *vtx = sim->msh->vtx.dat;

    double err = 0;
    int    num = 0;

    for (int h = 0; h < sim->msh->hxd.len; ++h) {
        struct hxd *hxd = &sim->msh->hxd.dat[h];

        int v0 = hxd->vtx[0];
        int v7 = hxd->vtx[7];

        double x0 = vtx[v0].dat[0];
        double y0 = vtx[v0].dat[1];
        double z0 = vtx[v0].dat[2];

        double x1 = vtx[v7].dat[0];
        double y1 = vtx[v7].dat[1];
        double z1 = vtx[v7].dat[2];

        int xn = (x1 - x0) / XS;
        int yn = (y1 - y0) / YS;
        int zn = (z1 - z0) / ZS;

        struct vec v;

        vec_new(&v, 3, x0, y0, z0);

        for (int i = 0; i < xn; ++i) {
            for (int j = 0; j < yn; ++j) {
                for (int k = 0; k < zn; ++k) {
                    double tgt = dat->tgt(&((struct sim_fun_ctx){.sim = sim}), &v);
                    double apx = sim->slv->apx(&((struct apx_fun_ctx){
                                                   .sim = sim,
                                                   .wgt = sim->slv->run.wgt[0],
                                                   .vtx = -1,
                                                   .qud = -1,
                                                   .hxd = h,
                                               }),
                        &v);

                    err += fabs(tgt - apx);
                    num += 1;

                    v.dat[2] += ZS;
                }

                v.dat[1] += YS;
            }

            v.dat[0] += XS;
        }
    }

    char name[256];

    fprintf(dat->fs, "%d & %.7e & %d & %.7e\n", sim->slv->run.ti, err / num, dat->non.itr, dat->non.err);

    if (sim->slv->ops.non.map && sim->slv->run.ti < sim->ops.tdd.num) {
        fclose(dat->fn);
        sprintf(name, "%s/non-%d.dat", sim->ops.exp.dir, sim->slv->run.ti + 1);
        dat->fn = fopen(name, "w+");
    }
}

void on_non(void *ctx, struct non_ops *ops)
{
    struct dat *dat = (struct dat *)ctx;

    dat->non.itr = ops->run.itr;
    dat->non.err = ops->run.err;

    double rlx = ops->run.rlx;

    fprintf(dat->fn, "%d & %.7e & %.3f\n", dat->non.itr, dat->non.err, rlx);
}

int pde(int argc, char **argv)
{
    if (argc < 3) {
        man();
        return -1;
    }

    struct sim sim;
    struct dat dat;

    if (sim_new(&sim)) {
        perror("fatal: initialization");
        return -1;
    }

    sim.slv->ops.non.mod = NON_NEW;
    sim.slv->ops.non.ops.dif = DIF_NUM;
    sim.slv->ops.non.ops.rlx = true;

    sprintf(sim.ops.usr.pfx, "u2");

    if (sim_imp_elm(&sim, argv[2])) {
        perror("fatal: import");
        return -1;
    }

    dat.tgt = dlsym(sim.ops.usr.hdl, "target");
    dat.non.itr = 0;
    dat.non.err = 0;

    sim.slv->itr.ctx = &dat;
    sim.slv->itr.run = on_slv;

    sim.slv->ops.non.ops.itr.ctx = &dat;
    sim.slv->ops.non.ops.itr.run = on_non;

    char name[256];

    sprintf(name, "%s/slv.dat", sim.ops.exp.dir);
    dat.fs = fopen(name, "w+");

    if (sim.slv->ops.non.map) {
        sprintf(name, "%s/non-0.dat", sim.ops.exp.dir);
        dat.fn = fopen(name, "w+");
    }

    if (sim_run(&sim)) {
        perror("fatal: execution");
        return -1;
    }

    fclose(dat.fs);

    if (sim.slv->ops.non.map) {
        fclose(dat.fn);
    }

    sim_cls(&sim);

    return 0;
}
