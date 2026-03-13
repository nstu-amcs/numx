#include <numx/vec/mtx.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int imtx_new(struct imtx *m, struct imtx_pps pps)
{
    assert(m);

    m->pps = pps;
    m->dat = malloc(sizeof(double *) * pps.r);

    if (!m->dat)
        return -1;

    for (int i = 0; i < pps.r; ++i) {
        m->dat[i] = malloc(sizeof(double) * pps.c);

        if (!m->dat[i]) {
            for (int j = 0; j < i; ++j)
                free(m->dat[j]);

            free(m->dat);

            return -1;
        }

        memset(m->dat[i], 0, sizeof(double) * pps.c);
    }

    return 0;
}

int imtx_cls(struct imtx *m)
{
    assert(m);

    for (int i = 0; i < m->pps.r; ++i)
        free(m->dat[i]);

    free(m->dat);

    return 0;
}

int imtx_vmul(struct imtx *m, struct vec *v, struct vec *r)
{
    assert(m);
    assert(v);
    assert(r);

    double **md = m->dat;
    double  *vd = v->dat;
    double  *rd = r->dat;

    // [[omp::directive(parallel for)]]
    for (int i = 0; i < m->pps.r; ++i) {
        double s = 0;

        // [[omp::directive(parallel for reduction(+ : s))]]
        for (int j = 0; j < m->pps.c; ++j)
            s += md[i][j] * vd[j];

        rd[i] = s;
    }

    return 0;
}

int imtx_mmul(struct imtx *a, struct imtx *b, struct imtx *r)
{
    assert(a);
    assert(b);
    assert(a);

    double **ad = a->dat;
    double **bd = b->dat;
    double **rd = r->dat;

    // [[omp::directive(parallel for)]]
    for (int i = 0; i < r->pps.r; ++i) {

        // [[omp::directive(parallel for)]]
        for (int j = 0; j < r->pps.c; ++j) {
            double s = 0;

            // [[omp::directive(parallel for reduction(+ : s))]]
            for (int e = 0; e < a->pps.c; ++e)
                s += ad[i][e] * bd[e][i];

            rd[i][j] = s;
        }
    }

    return 0;
}

int imtx_vdup(struct imtx *s, struct imtx *d)
{
    for (int i = 0; i < s->pps.r; ++i) {
        for (int j = 0; j < s->pps.c; ++j) {
            d->dat[i][j] = s->dat[i][j];
        }
    }

    return 0;
}
