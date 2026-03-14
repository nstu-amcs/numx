#include <math.h>
#include <stdio.h>

#include <numx/com/cmp.h>
#include <numx/msh/umsh.h>
#include <numx/non/apx.h>
#include <numx/pde/sim.h>

#define TOL 1e-10
#define MU0 (4 * M_PI * 1e-7)

#define BND_X_L -1
#define BND_X_R 1.07
#define BND_Y_B 0.00
#define BND_Y_T 1.05

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

static double lam_iron(void *c, struct vec *v)
{
    struct sim_fun_ctx *sim_fun_ctx = (struct sim_fun_ctx *)c;
    struct int_fun_ctx *int_fun_ctx = (struct int_fun_ctx *)sim_fun_ctx->ctx;
    struct sim         *sim = sim_fun_ctx->sim;

    struct apx_fun_ctx apx_ctx = {
        .sim = sim,
        .vtx = sim_fun_ctx->vtx,
        .seg = sim_fun_ctx->seg,
        .qud = sim_fun_ctx->qud,
        .var = 1,
    };

    double bx = sim->slv->apx.dif(&apx_ctx, v);
    apx_ctx.var = 0;
    double by = -sim->slv->apx.dif(&apx_ctx, v);
    double b = sqrt(bx * bx + by * by);
    double mu = int_cub_fun(int_fun_ctx, &(struct vec){.n = 1, .dat = &b});

    return 1.0 / (mu * MU0);
}

static void itr_cbk(void *, struct sim *sim)
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

        double az = sim->slv->apx.run(&ctx, &(struct vec){.n = 2, .dat = vtx->dat});
        ctx.var = 1;
        double bx = sim->slv->apx.dif(&ctx, &(struct vec){.n = 2, .dat = vtx->dat});
        ctx.var = 0;
        double by = -sim->slv->apx.dif(&ctx, &(struct vec){.n = 2, .dat = vtx->dat});

        printf("Az(%.2e, %.2e) = %.7e\n", vtx->dat[0], vtx->dat[1], az);
        printf("Bx(%.2e, %.2e) = %.7e\n", vtx->dat[0], vtx->dat[1], bx);
        printf("By(%.2e, %.2e) = %.7e\n", vtx->dat[0], vtx->dat[1], by);
        printf("|B|(%.2e, %.2e) = %.7e\n", vtx->dat[0], vtx->dat[1], sqrt(bx * bx + by * by));
    }

    struct vec bx;
    struct vec by;
    struct vec bm;

    vec_new(&bx, sim->msh->vtx.v2d.len);
    vec_new(&by, sim->msh->vtx.v2d.len);
    vec_new(&bm, sim->msh->vtx.v2d.len);

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        struct qud *qud = &sim->msh->qud.dat[qi];

        ctx.qud = qi;

        for (int vi = 0; vi < 4; ++vi) {
            int         vgi = qud->vtx[vi];
            struct v2d *vtx = &sim->msh->vtx.v2d.dat[vgi];

            ctx.vtx = vgi;
            ctx.var = 1;
            double _bx = sim->slv->apx.dif(&ctx, &(struct vec){.n = 2, .dat = vtx->dat});
            ctx.var = 0;
            double _by = -sim->slv->apx.dif(&ctx, &(struct vec){.n = 2, .dat = vtx->dat});

            bx.dat[vgi] = _bx;
            by.dat[vgi] = _by;
            bm.dat[vgi] = sqrt(_bx * _bx + _by * _by);
        }
    }

    sim->ops.exp.put_v(sim, "Bx", 1, &bx);
    sim->ops.exp.put_v(sim, "By", 2, &by);
    sim->ops.exp.put_v(sim, "Bm", 3, &bm);

    vec_cls(&bx);
    vec_cls(&by);
    vec_cls(&bm);
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

    if ((r = umsh_imp_tel(&msh, "msh", "msh"))) {
        goto end;
    }

    umsh_seg_srh(&msh, bnd_t, 0);
    umsh_seg_srh(&msh, bnd_l, 0);
    umsh_seg_srh(&msh, bnd_r, 0);
    umsh_seg_srh(&msh, bnd_b, 1);
    seg_cut_shr(&msh.seg);

    printf("[umsh][tel] seg: %d\n", msh.seg.len);

    if ((r = sim_imp_tel(&sim, "msh/msh.tel"))) {
        goto end;
    }

    // Create spline

    FILE *in = fopen("mu/mu.001", "r");
    int   n = 0;

    if (!in) {
        r = EIO;
        fclose(in);
        goto end;
    }

    fscanf(in, "%d", &n);

    struct imtx k;
    struct vec  x;

    if ((r = imtx_new(&k, (struct imtx_pps){.r = 4, .c = n}))) {
        fclose(in);
        goto end;
    }

    if ((r = vec_new(&x, n))) {
        fclose(in);
        goto end;
    }

    for (int i = 0; i < n; ++i) {
        fscanf(in, "%lf %lf", &k.dat[0][i], &x.dat[i]);
    }

    fclose(in);

    if ((r = apx_int_cub(&x, &k))) {
        goto end;
    }

    struct int_fun_ctx int_ctx = {
        .k = &k,
        .x = &x,
        .prv = -1,
    };

    sim.mat.dat[0].lam.type = VAL_FUN;
    sim.mat.dat[0].lam.as.fun.ctx = &int_ctx;
    sim.mat.dat[0].lam.as.fun.run = lam_iron;

    for (int i = 0; i < sim.mat.len; ++i) {
        if (sim.mat.dat[i].lam.type == VAL_NUM) {
            double mu = sim.mat.dat[i].lam.as.num;
            sim.mat.dat[i].lam.as.num = 1.0 / (mu * MU0);
        }
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
    strcpy(sim.ops.exp.dir, "out");
    strcpy(sim.ops.exp.pfx, "emf");
    strcpy(sim.ops.exp.sol, "Az");

    sim.slv->ops.iss.mod = ISS_GMR;
    sim.slv->ops.iss.ops.gmr.ops.err = 1e-7;
    sim.slv->ops.iss.ops.gmr.ops.itr.ctx = NULL;
    sim.slv->ops.iss.ops.gmr.ops.itr.run = NULL;
    sim.slv->ops.iss.ops.gmr.ops.max = 3000;

    sim.slv->itr_cbk.run = itr_cbk;
    sim.slv->ops.non.enable = true;

    if ((r = sim_run(&sim))) {
        goto end;
    }

end:
    umsh_cls(&msh);
    sim_cls(&sim);

    return r;
}
