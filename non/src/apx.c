#include <stdio.h>

#include <numx/non/apx.h>
#include <numx/vec/iss.h>

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

    printf("[apx][int-cub]|[iss-rlx] itr: %d\n", ops.ops.run.itr);
    printf("[apx][int-cub]|[iss-rlx] err: %.3e\n", ops.ops.run.err);

    for (int i = 0; i < n; ++i) {
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

end:
    mtx_cls(&mm);
    vec_cls(&gv);
    vec_cls(&hv);

    return r;
}
