#include <numx/non/fun.h>

double dif_tpm(void *ctx, mfun fun, struct vec *vtx, struct dif_ops *ops)
{
    assert(fun);
    assert(ops);

    int    var = ops->var;
    double hop = ops->hop != 0 ? ops->hop : NUMX_NON_DIF_HOP;

    if (ops->twk)
        ops->twk(ctx, hop, ops);
    else
        vtx->dat[var] += hop;

    double f1 = fun(ctx, vtx);

    if (ops->twk)
        ops->twk(ctx, -2 * hop, ops);
    else
        vtx->dat[var] -= 2 * hop;

    double f2 = fun(ctx, vtx);

    if (ops->twk)
        ops->twk(ctx, hop, ops);
    else
        vtx->dat[var] += hop;

    return (f1 - f2) / (2 * hop);
}
