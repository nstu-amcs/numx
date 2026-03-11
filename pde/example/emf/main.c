#include <math.h>
#include <stdio.h>

#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>
#include <numx/pde/sim.h>

#define TOL      1e-10
#define LAM_BASE (4 * M_PI * 1e-7)

#define BND_X_L -0.01
#define BND_X_R 0.08
#define BND_Y_B 0.00
#define BND_Y_T 0.06

static int bnd_t(struct umsh *msh, struct seg *s)
{
    struct v2d *a = &msh->vtx.v2d.dat[s->vtx[0]];
    struct v2d *b = &msh->vtx.v2d.dat[s->vtx[1]];

    const double ay = a->dat[1];
    const double by = b->dat[1];

    if (isclose(ay, BND_Y_T, TOL) && isclose(by, BND_Y_T, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_b(struct umsh *msh, struct seg *s)
{
    struct v2d *a = &msh->vtx.v2d.dat[s->vtx[0]];
    struct v2d *b = &msh->vtx.v2d.dat[s->vtx[1]];

    const double ay = a->dat[1];
    const double by = b->dat[1];

    if (isclose(ay, BND_Y_B, TOL) && isclose(by, BND_Y_B, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_l(struct umsh *msh, struct seg *s)
{
    struct v2d *a = &msh->vtx.v2d.dat[s->vtx[0]];
    struct v2d *b = &msh->vtx.v2d.dat[s->vtx[1]];

    const double ax = a->dat[0];
    const double bx = b->dat[0];

    if (isclose(ax, BND_X_L, TOL) && isclose(bx, BND_X_L, TOL)) {
        return 1;
    }

    return 0;
}

static int bnd_r(struct umsh *msh, struct seg *s)
{
    struct v2d *a = &msh->vtx.v2d.dat[s->vtx[0]];
    struct v2d *b = &msh->vtx.v2d.dat[s->vtx[1]];

    const double ax = a->dat[0];
    const double bx = b->dat[0];

    if (isclose(ax, BND_X_R, TOL) && isclose(bx, BND_X_R, TOL)) {
        return 1;
    }

    return 0;
}

static int cbk(void *, struct sim *sim)
{
    struct apx_fun_ctx ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
        .wgt = NULL,
    };

    struct v2d points[5] = {
        {.dat = {2.52e-2, 1.60e-3}},
        {.dat = {3.35e-2, 2.70e-3}},
        {.dat = {3.50e-2, 9.00e-4}},
        {.dat = {3.65e-2, 2.70e-3}},
        {.dat = {4.29e-2, 3.50e-3}},
    };

    for (int i = 0; i < 5; ++i) {
        struct v2d *vtx = &points[i];
        double      az = sim->slv->apx(&ctx, (union vtx_ptr){.v2d = vtx});

        printf("Az(%.2e, %.2e) = %.7e\n", vtx->dat[0], vtx->dat[1], az);
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

    if ((r = umsh_imp_tel(&msh, "pde/example/emf/msh", "msh"))) {
        goto end;
    }

    umsh_seg_srh(&msh, bnd_t, 0);
    umsh_seg_srh(&msh, bnd_l, 0);
    umsh_seg_srh(&msh, bnd_r, 0);
    umsh_seg_srh(&msh, bnd_b, 1);
    seg_cut_shr(&msh.seg);
    printf("Segments: %d\n", msh.seg.len);

    if ((r = sim_imp_tel(&sim, "pde/example/emf/msh/msh.tel"))) {
        goto end;
    }

    for (int i = 0; i < sim.mat.len; ++i) {
        double lam = sim.mat.dat[i].lam.as.num;
        sim.mat.dat[i].lam.as.num = 1.0 / (lam * LAM_BASE);
    }

    bnd_cut_dev(&sim.bnd, 2);
    sim.bnd.dat[0].cnd = 0;
    sim.bnd.dat[1].cnd = 1;

    cnd_bnd_cut_dev(&sim.cnd_bnd, 2);

    sim.cnd_bnd.dat[0].type = CND_BND_DIR;
    sim.cnd_bnd.dat[0].pps.dir.tgt.type = VAL_NUM;
    sim.cnd_bnd.dat[0].pps.dir.tgt.as.num = 0.0;

    sim.cnd_bnd.dat[1].type = CND_BND_NEU;
    sim.cnd_bnd.dat[1].pps.neu.tta.type = VAL_NUM;
    sim.cnd_bnd.dat[1].pps.neu.tta.as.num = 0.0;

    sim.slv = &slv.slv;
    sim.ops.exp.mod = SIM_EXP_CGNS;
    strcpy(sim.ops.exp.dir, "pde/example/emf/out");
    strcpy(sim.ops.exp.pfx, "emf");
    strcpy(sim.ops.exp.sol, "Magnetic Vector Potential");

    sim.slv->ops.iss.ops.bcg.ops.err = 1e-5;
    sim.slv->ops.iss.ops.bcg.ops.max = 2000;

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