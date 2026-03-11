#include <assert.h>
#include <math.h>
#include <numx/vec/mtx.h>
#include <stdlib.h>
#include <string.h>

int smtx_new(struct smtx *m, struct smtx_pps pps)
{
    assert(m);

    m->pps = pps;
    m->dr = malloc(sizeof(double) * pps.n);
    m->ia = malloc(sizeof(int) * (pps.n + 1));

    if (!m->dr || !m->ia) {
        free(m->dr);
        free(m->ia);

        return -1;
    }

    memset(m->dr, 0, sizeof(double) * pps.n);
    memset(m->ia, 0, sizeof(int) * (pps.n + 1));

    m->lr = NULL;
    m->ur = NULL;
    m->ja = NULL;

    if (pps.z > 0) {
        m->lr = malloc(sizeof(double) * pps.z);
        m->ur = malloc(sizeof(double) * pps.z);
        m->ja = malloc(sizeof(int) * pps.z);

        if (!m->lr || !m->ur || !m->ja) {
            free(m->dr);
            free(m->ia);
            free(m->lr);
            free(m->ur);
            free(m->ja);

            return -1;
        }

        memset(m->lr, 0, sizeof(double) * pps.z);
        memset(m->ur, 0, sizeof(double) * pps.z);
        memset(m->ja, 0, sizeof(int) * pps.z);
    }

    return 0;
}

int smtx_cls(struct smtx *m)
{
    assert(m);

    free(m->dr);
    free(m->lr);
    free(m->ur);
    free(m->ia);
    free(m->ja);

    return 0;
}

int smtx_inc(struct smtx *m, int i, int j, double v)
{
    assert(m);

    if (i < j) {
        int p = m->ia[j];

        while (m->ja[p] < i)
            ++p;

        if (m->ja[p] != i)
            return -1;

        m->ur[p] += v;

        return 0;
    }

    if (j < i) {
        int p = m->ia[i];

        while (m->ja[p] < j)
            ++p;

        if (m->ja[p] != j)
            return -1;

        m->lr[p] += v;

        return 0;
    }

    m->dr[i] += v;

    return 0;
}

int smtx_row_rst(struct smtx *m, int i)
{
    assert(m);

    m->dr[i] = 0.0;

    int p = m->ia[i];     // pointer to the first index in the target row
    int n = m->ia[i + 1]; // pointer to the first index in the next row

    while (p < n) {
        m->lr[p++] = 0.0;
    }

    for (int ui = 0; ui < m->pps.z; ++ui) {
        if (m->ja[ui] == i) {
            m->ur[ui] = 0.0;
        }
    }

    return 0;
}

int smtx_cmb(struct smtx *a, struct smtx *b, struct smtx *r, double k)
{
    assert(a);
    assert(b);
    assert(r);

    int n = a->pps.n;
    int z = a->pps.z;

    for (int i = 0; i < n; ++i)
        r->dr[i] = a->dr[i] + b->dr[i] * k;

    for (int i = 0; i < z; ++i) {
        r->lr[i] = a->lr[i] + b->lr[i] * k;
        r->ur[i] = a->ur[i] + b->ur[i] * k;
    }

    return 0;
}

int smtx_vdup(struct smtx *s, struct smtx *d)
{
    assert(s);
    assert(d);

    int n = s->pps.n;
    int z = s->pps.z;

    memcpy(d->ia, s->ia, sizeof(int) * (n + 1));
    memcpy(d->dr, s->dr, sizeof(double) * n);

    if (z > 0) {
        memcpy(d->ja, s->ja, sizeof(int) * z);
        memcpy(d->lr, s->lr, sizeof(double) * z);
        memcpy(d->ur, s->ur, sizeof(double) * z);
    }

    return 0;
}

int smtx_sdup(struct smtx *s, struct smtx *d)
{
    assert(s);
    assert(d);

    int n = s->pps.n;
    int z = s->pps.z;

    memcpy(d->ia, s->ia, sizeof(int) * (n + 1));

    if (z > 0)
        memcpy(d->ja, s->ja, sizeof(int) * z);

    return 0;
}

int smtx_rst(struct smtx *m)
{
    assert(m);

    int n = m->pps.n;
    int z = m->pps.z;

    memset(m->dr, 0, sizeof(double) * n);

    if (z > 0) {
        memset(m->lr, 0, sizeof(double) * z);
        memset(m->ur, 0, sizeof(double) * z);
    }

    return 0;
}

int smtx_ilu(struct smtx *m, struct smtx *r)
{
    assert(m);
    assert(r);

    int n = m->pps.n;
    int z = m->pps.z;

    int *mia = m->ia;
    int *mja = m->ja;

    int *ria = r->ia;
    int *rja = r->ja;

    memcpy(ria, mia, sizeof(int) * (n + 1));
    memcpy(rja, mja, sizeof(int) * z);

    double *mdr = m->dr;
    double *mlr = m->lr;
    double *mur = m->ur;

    double *rdr = r->dr;
    double *rlr = r->lr;
    double *rur = r->ur;

    for (int d = 0; d < n; ++d) {
        int lr0 = ria[d];
        int lr1 = ria[d + 1];

        double sd = 0;

        for (int lr = lr0; lr < lr1; ++lr) {
            int j = rja[lr];

            int ur0 = ria[j];
            int ur1 = ria[j + 1];

            double sl = 0;
            double su = 0;

            for (int lrr = lr0, urr = ur0; urr < ur1 && lrr < lr;) {
                int lj = rja[lrr];
                int ui = rja[urr];

                if (lj == ui) {
                    sl += rlr[lrr] * rur[urr];
                    su += rlr[urr] * rur[lrr];

                    lrr += 1;
                    urr += 1;

                    continue;
                }

                if (lj < ui)
                    lrr += 1;
                else
                    urr += 1;
            }

            rlr[lr] = (mlr[lr] - sl) / rdr[j];
            rur[lr] = (mur[lr] - su) / rdr[j];

            sd += rlr[lr] * rur[lr];
        }

        rdr[d] = sqrt(mdr[d] - sd);
    }

    return 0;
}

int smtx_dgl(struct smtx *m, struct smtx *r)
{
    assert(m);
    assert(r);

    int n = m->pps.n;

    double *mdr = m->dr;
    double *rdr = r->dr;

    for (int i = 0; i < n; ++i)
        rdr[i] = sqrt(mdr[i]);

    return 0;
}

int smtx_vmul(struct smtx *m, struct vec *x, struct vec *f)
{
    assert(m);
    assert(x);
    assert(f);

    int n = m->pps.n;

    int *ia = m->ia;
    int *ja = m->ja;

    double *dr = m->dr;
    double *lr = m->lr;
    double *ur = m->ur;

    double *xv = x->dat;
    double *fv = f->dat;

    for (int i = 0; i < n; ++i)
        fv[i] = xv[i] * dr[i];

    for (int i = 0; i < n; ++i) {
        int ar0 = ia[i];
        int ar1 = ia[i + 1];

        for (int ar = ar0; ar < ar1; ++ar) {
            int j = ja[ar];

            fv[i] += xv[j] * lr[ar];
            fv[j] += xv[i] * ur[ar];
        }
    }

    return 0;
}
