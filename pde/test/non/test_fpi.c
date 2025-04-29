#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <numx/pde/sim.h>

static const char *usr[] = {"u1", "u2", "u3", "u4", "u5", NULL};

struct ctx
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

static void on_slv(void *ctx, struct sim *sim);
static void on_non(void *ctx, struct non_ops *ops);

int main()
{
    struct sim sim;
    struct ctx ctx = {
        .fs = fopen("pde/test/non/case/res/fpi.slv", "w+"),
        .fn = NULL,
        .non = {
                .itr = 0,
                .err = 0,
                }
    };

    char fn[256];

    for (int i = 0; usr[i] != NULL; ++i) {
        sim_new(&sim);

        strcpy(sim.ops.usr.pfx, usr[i]);
        strcpy(sim.ops.exp.pfx, usr[i]);

        sim_imp_elm(&sim, "pde/test/non/case/elm/case.sif");

        sim.slv->ops.non.ops.itr.ctx = &ctx;
        sim.slv->ops.non.ops.itr.run = on_non;
        sim.slv->itr.ctx = &ctx;
        sim.slv->itr.run = on_slv;

        ctx.tgt = (double (*)(struct sim_fun_ctx *, struct vec *))dlsym(sim.ops.usr.hdl, "target");

        sim.slv->ops.non.mod = NON_FPI;
        sim.slv->ops.non.ops.rlx = false;

        sprintf(fn, "%s/fpi.%s.non", sim.ops.exp.dir, usr[i]);
        ctx.fn = fopen(fn, "w+");

        sim_run(&sim);
        sim_cls(&sim);
    }

    return 0;
}

static void on_slv(void *c, struct sim *sim)
{
    struct ctx *ctx = (struct ctx *)c;
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

    struct vec v;
    vec_new(&v, 3);

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
                    double tgt = ctx->tgt(&((struct sim_fun_ctx){.sim = sim}), &v);
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

    fprintf(ctx->fs, "%d & %.7e & %.7e\n", ctx->non.itr, err / num, ctx->non.err);
}

static void on_non(void *c, struct non_ops *ops)
{
    struct ctx *ctx = (struct ctx *)c;

    ctx->non.itr = ops->run.itr;
    ctx->non.err = ops->run.err;
    double rlx = ops->run.rlx;

    fprintf(ctx->fn, "%d & %.7e & %.3f\n", ctx->non.itr, ctx->non.err, rlx);
}
