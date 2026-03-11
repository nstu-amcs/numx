#include <stdio.h>

#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>
#include <numx/pde/sim.h>

#define TOL 1e-7

static double u(void *, struct vec *v)
{
    const double x = v->dat[0];
    const double y = v->dat[1];

    return x * x + y + 5.0;
}

static double f(void *, struct vec *v)
{
    const double x = v->dat[0];
    const double y = v->dat[1];

    return y;
}

static double tta(void *, struct vec *v)
{
    const double x = v->dat[0];
    const double y = v->dat[1];

    return -4.0 * x;
}

static int bnd_t(vtx_ptr a, vtx_ptr b)
{
    const double ay = a.v2d->dat[1];
    const double by = b.v2d->dat[1];

    if (isclose(ay, 1.0, TOL) && isclose(by, 1.0, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_b(vtx_ptr a, vtx_ptr b)
{
    const double ay = a.v2d->dat[1];
    const double by = b.v2d->dat[1];

    if (isclose(ay, 0.0, TOL) && isclose(by, 0.0, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_l(vtx_ptr a, vtx_ptr b)
{
    const double ax = a.v2d->dat[0];
    const double bx = b.v2d->dat[0];

    if (isclose(ax, 0.0, TOL) && isclose(bx, 0.0, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_r(vtx_ptr a, vtx_ptr b)
{
    const double ax = a.v2d->dat[0];
    const double bx = b.v2d->dat[0];

    if (isclose(ax, 2.0, TOL) && isclose(bx, 2.0, TOL)) {
        return 1;
    }

    return 0;
}

static int cbk(void *, struct sim *sim)
{
    for (int i = 0; i < sim->msh->vtx.v2d.len; ++i) {
        printf("%.3lf ", sim->slv->run.wgt[0]->dat[i]);
    }

    printf("\n");

    for (int i = 0; i < sim->msh->vtx.v2d.len; ++i) {
        struct v2d *vtx = &sim->msh->vtx.v2d.dat[i];
        printf("%.3lf ", u(NULL, &(struct vec){.n = 2, .dat = vtx->dat}));
    }

    printf("\n");

    return 0;
}

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    int r = 0;

    struct umsh msh;
    struct fem  slv;
    struct sim  sim;

    if ((r = umsh_new(&msh))) {
        goto end;
    }

    if ((r = fem_new(&slv))) {
        goto end;
    }

    if ((r = sim_new(&sim))) {
        goto end;
    }

    sim.msh = &msh;

    if ((r = umsh_imp_tel(&msh, "pde/example/p2-u/msh", "msh"))) {
        goto end;
    }

    umsh_seg_srh(&msh, bnd_b, 0);
    umsh_seg_srh(&msh, bnd_r, 0);
    umsh_seg_srh(&msh, bnd_t, 0);
    umsh_seg_srh(&msh, bnd_l, 1);
    seg_cut_shr(&msh.seg);
    printf("Segments: %d\n", msh.seg.len);

    if ((r = sim_imp_tel(&sim, "pde/example/p2-u/msh/msh.tel"))) {
        goto end;
    }

    bnd_cut_dev(&sim.bnd, 2);
    sim.bnd.dat[0].cnd = 0;
    sim.bnd.dat[1].cnd = 1;

    cnd_bnd_cut_dev(&sim.cnd_bnd, 2);

    sim.cnd_bnd.dat[0].type = CND_BND_DIR;
    sim.cnd_bnd.dat[0].pps.dir.tgt.type = VAL_FUN;
    sim.cnd_bnd.dat[0].pps.dir.tgt.as.fun = u;

    // sim.cnd_bnd.dat[1].type = CND_BND_NEU;
    // sim.cnd_bnd.dat[1].pps.neu.tta.type = VAL_FUN;
    // sim.cnd_bnd.dat[1].pps.neu.tta.as.fun = tta;

    sim.cnd_bnd.dat[1].type = CND_BND_NEU;
    sim.cnd_bnd.dat[1].pps.neu.tta.type = VAL_NUM;
    sim.cnd_bnd.dat[1].pps.neu.tta.as.num = 0.0;

    sim.slv = &slv.slv;
    sim.ops.exp.mod = SIM_EXP_CGNS;
    strcpy(sim.ops.exp.dir, "pde/example/p2-u/out");
    strcpy(sim.ops.exp.pfx, "temp");
    strcpy(sim.ops.exp.sol, "Temperature");

    sim.slv->itr_cbk.run = cbk;

    if ((r = sim_run(&sim))) {
        goto end;
    }

    printf("Error: %.7lf\n", sim.slv->ops.iss.ops.bcg.ops.run.err);
    printf("Iterations: %d\n", sim.slv->ops.iss.ops.bcg.ops.run.itr);

end:
    umsh_cls(&msh);
    sim_cls(&sim);

    return r;
}