#include <stdio.h>

#include <numx/non/apx.h>
#include <numx/vec/iss.h>
#include <string.h>

int apx_int_cub(struct vec *xv, struct imtx *km)
{
    if (!xv || !km) {
        return -EINVAL;
    }

    int n = xv->n - 1; // number of subintervals
    int r = 0;

    struct dmtx mm;
    struct vec  gv;
    struct vec  hv;

    if ((r = dmtx_new(&mm, (struct dmtx_pps){.n = n + 1, .d = 3}))) {
        goto end;
    }

    mm.la[0] = 0;
    mm.la[1] = -1;
    mm.la[2] = 1;

    if ((r = vec_new(&gv, n + 1))) {
        goto end;
    }

    if ((r = vec_new(&hv, n))) {
        goto end;
    }

    double *g = gv.dat;
    double *h = hv.dat;
    double *x = xv->dat;

    double *a = km->dat[0];
    double *b = km->dat[1];
    double *c = km->dat[2];
    double *d = km->dat[3];

    for (int i = 0; i < n; ++i)
        h[i] = x[i + 1] - x[i];

    for (int i = 1; i < n; ++i)
        g[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

    g[0] = 0;
    g[n] = 0;

    mm.ad[0][0] = 1;
    mm.ad[n][0] = 1;

    for (int i = 1; i < n; ++i) {
        mm.ad[i][0] = 2 * (h[i - 1] + h[i]);
        mm.ad[i][1] = h[i - 1];
        mm.ad[i][2] = h[i];
    }

    struct vec cv = {
        .n = n + 1,
        .dat = c,
    };

    struct iss_rlx_ops ops = {
        .ops =
            {
                  .err = 1e-10,
                  .max = 1000,
                  .itr =
                    {
                        .ctx = NULL,
                        .run = NULL,
                    }, },
        .rlx = 1,
    };

    if ((r = diss_rlx_slv(&mm, &cv, &gv, &ops))) {
        goto end;
    }

    for (int i = 0; i < n; ++i) {
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    printf("[apx][int-cub] ok\n");

end:
    mtx_cls(&mm);
    vec_cls(&gv);
    vec_cls(&hv);

    if (r) {
        printf("[apx][int-cub] err: %s\n", strerror(-r));
    }

    return r;
}

double int_cub_fun(void *ctx, struct vec *x)
{
    struct int_fun_ctx *int_fun_ctx = (struct int_fun_ctx *)ctx;
    struct vec         *ip = int_fun_ctx->x;
    struct imtx        *km = int_fun_ctx->k;

    double v = x->dat[0];
    int    i = int_fun_ctx->prv;
    int    n = int_fun_ctx->x->n;

    if (i == -1) {
        i = 0;
    }

    double x0 = ip->dat[i];
    double x1 = ip->dat[i + 1];

    while (v < x0) {
        if (i == 0) {
            // Hit the first interval, perform linear approximation: y = ax + b.

            double y0 = km->dat[0][i];
            double y1 = km->dat[0][i + 1];
            double a = (y1 - y0) / (x1 - x0);
            double b = y0 - a * x0;

            int_fun_ctx->prv = i;

            return a * v + b;
        }

        i -= 1;
        x0 = ip->dat[i];
        x1 = ip->dat[i + 1];
    }

    while (v > x1) {
        if (i == n - 2) {
            // Hit the last interval, perform linear approximation: y = ax + b.

            double y0 = km->dat[0][i];
            double y1 = km->dat[0][i + 1];
            double a = (y1 - y0) / (x1 - x0);
            double b = y0 - a * x0;

            int_fun_ctx->prv = i;

            return a * v + b;
        }

        i += 1;
        x0 = ip->dat[i];
        x1 = ip->dat[i + 1];
    }

    double a = km->dat[0][i];
    double b = km->dat[1][i];
    double c = km->dat[2][i];
    double d = km->dat[3][i];
    double r = v - x0;

    int_fun_ctx->prv = i;

    return a + b * r + c * r * r + d * r * r * r;
}

double int_cub_dif(void *ctx, mfun fun, struct vec *vtx, struct dif_ops *ops)
{
    (void)fun;
    (void)ops;

    struct int_fun_ctx *int_fun_ctx = (struct int_fun_ctx *)ctx;
    struct vec         *ip = int_fun_ctx->x;
    struct imtx        *km = int_fun_ctx->k;

    double v = vtx->dat[0];
    int    i = int_fun_ctx->prv;
    int    n = int_fun_ctx->x->n;

    if (i == -1) {
        i = 0;
    }

    double x0 = ip->dat[i];
    double x1 = ip->dat[i + 1];

    while (v < x0) {
        if (i == 0) {
            // Hit the first interval, perform linear approximation: y = ax + b.

            double y0 = km->dat[0][i];
            double y1 = km->dat[0][i + 1];
            double a = (y1 - y0) / (x1 - x0);

            int_fun_ctx->prv = i;

            return a;
        }

        i -= 1;
        x0 = ip->dat[i];
        x1 = ip->dat[i + 1];
    }

    while (v > x1) {
        if (i == n - 2) {
            // Hit the last interval, perform linear approximation: y = ax + b.

            double y0 = km->dat[0][i];
            double y1 = km->dat[0][i + 1];
            double a = (y1 - y0) / (x1 - x0);

            int_fun_ctx->prv = i;

            return a;
        }

        i += 1;
        x0 = ip->dat[i];
        x1 = ip->dat[i + 1];
    }

    double b = km->dat[1][i];
    double c = km->dat[2][i];
    double d = km->dat[3][i];

    int_fun_ctx->prv = i;

    return b + c * 2 * (v - x0) + d * 3 * (v * v - 2 * v * x0 + x0 * x0);
}
