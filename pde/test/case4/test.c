#include <assert.h>
#include <math.h>
#include <numx/pde/sim.h>
#include <stdarg.h>
#include <stdio.h>

struct non_itr
{
    int    itr;
    double err;
};

static double target(struct vtx *v)
{
    return v->x + v->y + v->z;
}

static void sim_cbk(void *ctx, int n, ...)
{
    va_list arg;
    va_start(arg, n);

    double     *err = (double *)ctx;
    struct sim *sim = va_arg(arg, struct sim *);
    int         tm = va_arg(arg, int);
    struct vec *wgt = va_arg(arg, struct vec *);

    va_end(arg);

    (void)tm;

    sim_err(sim, wgt, target, err);
}

static void non_cbk(void *ctx, int n, ...)
{
    va_list arg;
    va_start(arg, n);

    struct non_itr *hld = (struct non_itr *)ctx;
    struct sim     *sim = va_arg(arg, struct sim *);
    int             tm = va_arg(arg, int);
    int             itr = va_arg(arg, int);
    double          err = va_arg(arg, double);

    va_end(arg);

    (void)sim;
    (void)tm;

    hld->itr = itr;
    hld->err = err;
}

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    struct sim sim;

    if (sim_new(&sim))
        return -1;

    if (sim_imp_elm(&sim, "case.sif"))
        return -1;

    struct non_itr non;
    double         err = 0;

    ((struct fem *)sim.slv)->ops.non.itr.ctx = &non;
    ((struct fem *)sim.slv)->ops.non.itr.run = non_cbk;

    ((struct fem *)sim.slv)->itr.ctx = &err;
    ((struct fem *)sim.slv)->itr.run = sim_cbk;

    if (sim_run(&sim))
        return -1;

    sim_cls(&sim);
    printf("Nonlinear: err = %lf, itr = %d\n", fabs(non.err), non.itr);
    printf("Simulation: err = %lf\n", fabs(err));

    return 0;
}
