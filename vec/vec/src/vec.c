#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <numx/vec/vec.h>

int vec_new_ini(struct vec *v, int n, ...)
{
    assert(v);
    assert(n > 0);

    v->n = n;
    v->dat = malloc(sizeof(double) * n);

    if (!v->dat)
        return -1;

    memset(v->dat, 0, sizeof(double) * n);

    va_list arg;
    va_start(arg, n);

    for (int i = 0; i < n; ++i) {
        double a = va_arg(arg, double);

        if (a == DBL_MAX)
            break;

        v->dat[i] = a;
    }

    va_end(arg);

    return 0;
}

int vec_cls(struct vec *v)
{
    assert(v);

    free(v->dat);

    v->n = 0;
    v->dat = NULL;

    return 0;
}

int vec_dst(struct vec *a, struct vec *b, double *r)
{
    assert(a);
    assert(b);
    assert(r);

    double s = 0;
    double p = 0;

    int n = a->n;

    double *ad = a->dat;
    double *bd = b->dat;

#ifdef OMP
    [[omp::directive(parallel for reduction(+:s))]]
#endif
    for (int i = 0; i < n; ++i) {
        p = bd[i] - ad[i];
        s += p * p;
    }

    *r = sqrt(s);

    return 0;
}

int vec_cmb(struct vec *a, struct vec *b, struct vec *r, double k)
{
    assert(a);
    assert(b);
    assert(r);

    assert(a->n == b->n);
    assert(b->n == r->n);

    int n = a->n;

    double *ad = a->dat;
    double *bd = b->dat;
    double *rd = r->dat;

#ifdef OMP
    [[omp::directive(parallel for)]]
#endif
    for (int i = 0; i < n; ++i)
        rd[i] = ad[i] + k * bd[i];

    return 0;
}

int vec_dot(struct vec *a, struct vec *b, double *r)
{
    assert(a);
    assert(b);
    assert(r);

    int dim = a->n;

    double *ad = a->dat;
    double *bd = b->dat;

    double s = 0;

#ifdef OMP
    [[omp::directive(parallel for reduction(+:s))]]
#endif
    for (int i = 0; i < dim; ++i)
        s += ad[i] * bd[i];

    *r = s;

    return 0;
}

int vec_nrm(struct vec *v, double *r)
{
    assert(v);
    assert(r);

    vec_dot(v, v, r);
    *r = sqrt(*r);

    return 0;
}

int vec_mul(struct vec *v, struct vec *r, double m)
{
    int n = v->n;

    double *vd = v->dat;
    double *rd = r->dat;

#ifdef OMP
    [[omp::directive(parallel for)]]
#endif
    for (int i = 0; i < n; ++i)
        rd[i] = vd[i] * m;

    return 0;
}

int vec_rot(
    struct vec *restrict v, struct vec *restrict r, int o, double c, double s)
{
    assert(v);
    assert(r);

    int n = v->n;
    int i = 0;

    for (i = 0; i < o; ++i)
        r->dat[i] = v->dat[i];

    r->dat[i] = v->dat[i] * c + v->dat[i + 1] * s;
    i += 1;
    r->dat[i] = -v->dat[i - 1] * s + v->dat[i] * c;

    for (i = i + 1; i < n; ++i)
        r->dat[i] = v->dat[i];

    return 0;
}

int vec_dup(struct vec *s, struct vec *d)
{
    assert(s);
    assert(d);

    memcpy(d->dat, s->dat, sizeof(double) * s->n);

    return 0;
}

int vec_swp(struct vec *a, struct vec *b)
{
    assert(a);
    assert(b);

    double *t = a->dat;
    a->dat = b->dat;
    b->dat = t;

    return 0;
}

int vec_rst(struct vec *v)
{
    assert(v);

    memset(v->dat, 0, sizeof(double) * v->n);

    return 0;
}

cut_gen(vec_cut, vec, PUB);
log_gen(vec_log, vec_rec, vec, PUB);
que_gen(vec_que, vec_cut, vec, PUB);
