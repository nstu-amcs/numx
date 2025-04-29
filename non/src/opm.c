#include <numx/non/opm.h>

double opm_loc_bis(void *ctx, double (*fun)(void *, struct vec *), struct opm_ops *ops)
{
    assert(fun);
    assert(ops);

    struct vec *vtx = ops->vtx;
    int         var = ops->var;
    double      eps = ops->eps != 0 ? ops->eps : NUMX_NON_OPM_EPS;
    double      del = eps / 2;

    double a = ops->beg;
    double b = ops->end;

    while (b - a > eps) {
        double x1 = (a + b - del) / 2;
        double x2 = (a + b + del) / 2;

        if (ops->twk)
            ops->twk(ctx, x1, ops);
        else
            vtx->dat[var] = x1;

        double f1 = fun(ctx, vtx);

        if (ops->twk)
            ops->twk(ctx, x2, ops);
        else
            vtx->dat[var] = x2;

        double f2 = fun(ctx, vtx);

        if (f1 < f2)
            b = x2;
        else
            a = x1;
    }

    return (a + b) / 2;
}
