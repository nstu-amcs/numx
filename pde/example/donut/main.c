#include <stdio.h>

#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>
#include <numx/pde/sim.h>

#define TOL 1e-7

static double tgt(void *, struct vec *v)
{
    // const double x = v->dat[0];
    const double y = v->dat[1];

    return y;
}

static int bnd_t(vtx_ptr a, vtx_ptr b)
{
    const double ay = a.v2d->dat[1];
    const double by = b.v2d->dat[1];

    if (isclose(ay, 3.0, TOL) && isclose(by, 3.0, TOL)) {
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

    if (isclose(ax, 3.0, TOL) && isclose(bx, 3.0, TOL)) {
        return 1;
    }

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

    if ((r = umsh_imp_tel(&msh, "pde/example/donut/msh", "msh"))) {
        goto end;
    }

    umsh_seg_srh(&msh, bnd_b, 0);
    umsh_seg_srh(&msh, bnd_l, 0);
    umsh_seg_srh(&msh, bnd_r, 0);
    umsh_seg_srh(&msh, bnd_t, 0);
    seg_cut_shr(&msh.seg);

    printf("Segmets: %d\n", msh.seg.len);

    if ((r = sim_imp_tel(&sim, "pde/example/donut/msh/msh.tel"))) {
        goto end;
    }

    bnd_cut_dev(&sim.bnd, 1);
    sim.bnd.dat[0].cnd = 0;

    cnd_bnd_cut_dev(&sim.cnd_bnd, 1);

    sim.cnd_bnd.dat[0].type = CND_BND_DIR;
    sim.cnd_bnd.dat[0].pps.dir.tgt.type = VAL_NUM;
    sim.cnd_bnd.dat[0].pps.dir.tgt.as.num = 0;

    sim.slv = &slv.slv;
    sim.ops.exp.mod = SIM_EXP_CGNS;
    strcpy(sim.ops.exp.dir, "pde/example/donut/out");
    strcpy(sim.ops.exp.pfx, "temp");
    strcpy(sim.ops.exp.sol, "Temperature");

    if ((r = sim_run(&sim))) {
        goto end;
    }

end:
    umsh_cls(&msh);
    sim_cls(&sim);

    return r;
}