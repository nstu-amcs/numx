#include <dlfcn.h>
#include <math.h>
#include <munit.h>
#include <stdio.h>
#include <stdlib.h>

#include <numx/pde/sim.h>

MunitSuite unit_pde_ell_suite();
MunitSuite unit_pde_pbc_suite();
MunitSuite unit_pde_hmc_suite();

struct ctx
{
    mfun  tgt;
    FILE *out;
};

void export(void *und, struct sim *sim)
{
    struct ctx *ctx = (struct ctx *)und;
    struct vec *vtx = sim->msh->vtx.dat;

    double err = 0;
    int    num = 0;

    struct apx_fun_ctx apx_ctx = {
        .sim = sim,
        .wgt = sim->slv->run.wgt[0],
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    double dat[3] = {0, 0, 0};

    struct vec v = {
        .n = 3,
        .dat = dat,
    };

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
    };

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

        double xs = (x1 - x0) / 2;
        double ys = (y1 - y0) / 2;
        double zs = (z1 - z0) / 2;

        apx_ctx.hxd = h;

        v.dat[0] = x0;
        v.dat[1] = y0;
        v.dat[2] = z0;

        for (int i = 0; i <= 2; ++i) {
            for (int j = 0; j <= 2; ++j) {
                for (int k = 0; k <= 2; ++k) {
                    double tgt = ctx->tgt(&fun_ctx, &v);
                    double apx = sim->slv->apx(&apx_ctx, &v);

                    err += fabs(tgt - apx);
                    num += 1;

                    v.dat[2] += zs;
                }

                v.dat[1] += ys;
            }

            v.dat[0] += xs;
        }
    }

    fprintf(ctx->out, "%d & %.7e\n", sim->slv->run.ti, err / num);
}

MunitResult test(const MunitParameter pps[], void *dir)
{
    char *slv = "fem";
    char *msh = "m1";
    char *frq = "f1";
    char *usr = "u1";
    char *bdf = "2";
    char *ini = "1";

    for (int i = 0; i < 6; ++i) {
        if (!strcmp("slv", pps[i].name)) {
            slv = pps[i].value;
            continue;
        }

        if (!strcmp("msh", pps[i].name)) {
            msh = pps[i].value;
            continue;
        }

        if (!strcmp("frq", pps[i].name)) {
            frq = pps[i].value;
            continue;
        }

        if (!strcmp("usr", pps[i].name)) {
            usr = pps[i].value;
            continue;
        }

        if (!strcmp("bdf", pps[i].name)) {
            bdf = pps[i].value;
            continue;
        }

        if (!strcmp("ini", pps[i].name)) {
            ini = pps[i].value;
            continue;
        }
    }

    char path[256];
    char exec[256];

    struct sim sim;

    if (sim_new(&sim))
        return MUNIT_FAIL;

    sim.ops.exp.mod = SIM_EXP_GNS;
    sim.ops.exp.ini = sim_exp_gns_ini;
    sim.ops.exp.put = sim_exp_gns_put;

    sprintf(sim.ops.msh.dir, "pde/test/msh/%s/%s/elm", msh, frq);
    sprintf(sim.ops.msh.pfx, "mesh");
    sprintf(sim.ops.usr.dir, "%s", (char *)dir);
    sprintf(sim.ops.usr.pfx, "%s", usr);
    sprintf(sim.ops.exp.dir, "%s/res", (char *)dir);
    sprintf(sim.ops.exp.pfx, "%s.%s.b%s", frq, usr, bdf);
    sprintf(exec,
        "gcc -o /tmp/numx_usr.so -I/usr/share/include -shared -fPIC %s/%s.c",
        sim.ops.usr.dir, sim.ops.usr.pfx);

    system(exec);

    if ((sim.ops.usr.hdl = dlopen("/tmp/numx_usr.so", RTLD_NOW)) == NULL)
        return MUNIT_FAIL;

    sim.msh = malloc(sizeof(struct msh));

    if (sim.msh == NULL || msh_new(sim.msh))
        return MUNIT_FAIL;

    if (msh_imp_grd(sim.msh, sim.ops.msh.dir, sim.ops.msh.pfx))
        return MUNIT_FAIL;

    if (!strcmp("fem", slv)) {
        struct fem *fem = malloc(sizeof(struct fem));

        if (fem == NULL || fem_new(fem))
            return MUNIT_FAIL;

        sim.slv = &fem->slv;
    }

    sprintf(path, "%s/elm.sif", (char *)dir);

    if (sim_imp_elm(&sim, path))
        return MUNIT_FAIL;

    sprintf(path, "%s/res/%s.%s.b%s.slv", (char *)dir, frq, usr, bdf);

    struct ctx ctx = {
        .tgt = dlsym(sim.ops.usr.hdl, "target"),
        .out = fopen(path, "w+"),
    };

    sim.slv->itr.ctx = &ctx;
    sim.slv->itr.run = export;
    sim.slv->ops.tdd = atoi(bdf);
    sim.ops.tdd.ini = atoi(ini);

    if (sim_run(&sim))
        return MUNIT_FAIL;

    fclose(ctx.out);
    sim_cls(&sim);

    return MUNIT_OK;
}

int main(int argc, char **argv)
{
    MunitSuite suites[] = {
        unit_pde_ell_suite(),
        unit_pde_pbc_suite(),
        unit_pde_hmc_suite(),
    };

    const MunitSuite suite = {
        "/numx/pde",
        NULL,                    // tests
        suites,                  // suites
        1,                       // iterations
        MUNIT_SUITE_OPTION_NONE, // options
    };

    return munit_suite_main(&suite, NULL, argc, argv);
}
