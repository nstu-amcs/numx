#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>
#include <numx/pde/sim.h>

#define TOL 1e-7

int bnd_t(vtx_ptr a, vtx_ptr b)
{
    double ax = a.v2d->dat[0];
    double ay = a.v2d->dat[1];
    double bx = b.v2d->dat[0];
    double by = b.v2d->dat[1];

    if (isclose(ay, 0.1, TOL) && isclose(by, 0.1, TOL)) {
        if (ax > bx) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
}

int bnd_b(vtx_ptr a, vtx_ptr b)
{
    double ax = a.v2d->dat[0];
    double ay = a.v2d->dat[1];
    double bx = b.v2d->dat[0];
    double by = b.v2d->dat[1];

    if (isclose(ay, 0.0, TOL) && isclose(by, 0.0, TOL)) {
        if (ax < bx) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
}

int bnd_l(vtx_ptr a, vtx_ptr b)
{
    double ax = a.v2d->dat[0];
    double ay = a.v2d->dat[1];
    double bx = b.v2d->dat[0];
    double by = b.v2d->dat[1];

    if (isclose(ax, -0.1, TOL) && isclose(bx, -0.1, TOL)) {
        if (ay > by) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
}

int bnd_r(vtx_ptr a, vtx_ptr b)
{
    double ax = a.v2d->dat[0];
    double ay = a.v2d->dat[1];
    double bx = b.v2d->dat[0];
    double by = b.v2d->dat[1];

    if (isclose(ax, 0.17, TOL) && isclose(bx, 0.17, TOL)) {
        if (ay < by) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
}

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    int r = 0;

    struct umsh msh;
    struct sim  sim;

    if (umsh_new(&msh)) {
        r = 1;
        goto end;
    }

    if (umsh_imp_tel(&msh, "pde/example/emf/msh", "msh")) {
        r = 1;
        goto end;
    }

    umsh_seg_srh(&msh, bnd_t, 0);
    umsh_seg_srh(&msh, bnd_l, 0);
    umsh_seg_srh(&msh, bnd_r, 0);
    umsh_seg_srh(&msh, bnd_b, 1);

    seg_cut_shr(&msh.seg);

    if (sim_new(&sim)) {
        r = 1;
        goto end;
    }

    if (sim_imp_tel(&sim, "pde/example/emf/msh/msh.tel")) {
        r = 1;
        goto end;
    }

    sim.msh = &msh;

    bnd_cut_dev(&sim.bnd, 2);
    sim.bnd.dat[0].cnd = 0;
    sim.bnd.dat[1].cnd = 1;

    cnd_bnd_cut_dev(&sim.cnd_bnd, 2);

    sim.cnd_bnd.dat[0].type = CND_BND_NEU;
    sim.cnd_bnd.dat[0].pps.neu.tta.type = VAL_NUM;
    sim.cnd_bnd.dat[0].pps.neu.tta.as.num = 0;

    sim.cnd_bnd.dat[1].type = CND_BND_DIR;
    sim.cnd_bnd.dat[1].pps.dir.tgt.type = VAL_NUM;
    sim.cnd_bnd.dat[1].pps.dir.tgt.as.num = 0;

    struct fem slv;

    if (fem_new(&slv)) {
        r = 1;
        goto end;
    }

    sim.slv = &slv.slv;
    strcpy(sim.ops.exp.dir, "pde/example/emf/out");
    strcpy(sim.ops.exp.pfx, "emf");
    sim.ops.exp.mod = SIM_EXP_GNS;

    if (sim_run(&sim)) {
        r = 1;
        goto end;
    }

end:
    umsh_cls(&msh);
    sim_cls(&sim);

    return r;
}