#include <numx/non/dif.h>

double pdif(void *ctx, double (*fun)(void *, struct vec *), int var, double hop, struct vec *vtx)
{
    assert(fun);
    assert(vtx);

    vtx->dat[var] += hop;
    double f1 = fun(ctx, vtx);

    vtx->dat[var] -= 2 * hop;
    double f2 = fun(ctx, vtx);

    vtx->dat[var] += hop;

    return (f1 - f2) / (2 * hop);
}
