#include <errno.h>
#include <math.h>
#include <string.h>

#include <numx/vec/iss.h>

static int siss_bcg_unc_slv(
    struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
{
    int n = m->pps.n;
    int c = 0;

    struct vec r;
    struct vec z;
    struct vec p;
    struct vec s;
    struct vec h;

    if (vec_new(&r, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&z, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&p, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&s, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&h, n)) {
        c = -1;
        goto end;
    }

    double alf = 0;
    double bet = 0;
    double omg = 0;
    double tmp = 0;
    double nrm = 0;

    mtx_vmul(m, x, &r);
    vec_cmb(f, &r, &r, -1);

    vec_dup(&r, &z);
    vec_dup(&r, &p);

    for (int k = 1; k <= o->ops.max; ++k) {
        mtx_vmul(m, &p, &h);

        vec_dot(&r, &z, &tmp);
        vec_dot(&h, &z, &bet);

        alf = tmp / bet;

        vec_cmb(&r, &h, &s, -alf);

        mtx_vmul(m, &s, &r);
        vec_dot(&r, &s, &omg);
        vec_dot(&r, &r, &bet);

        omg = omg / bet;

        vec_cmb(x, &p, x, alf);
        vec_cmb(x, &s, x, omg);
        vec_cmb(&s, &r, &r, -omg);
        vec_nrm(&r, &nrm);

        o->ops.run.itr = k;
        o->ops.run.err = nrm;

        if (o->ops.itr.run)
            o->ops.itr.run(o->ops.itr.ctx, &o->ops);

        if (nrm < o->ops.err)
            break;

        vec_dot(&r, &z, &bet);

        bet = (bet / tmp) * (alf / omg);

        vec_cmb(&p, &h, &p, -omg);
        vec_cmb(&r, &p, &p, bet);
    }

end:
    vec_cls(&r);
    vec_cls(&z);
    vec_cls(&p);
    vec_cls(&s);
    vec_cls(&h);

    return c;
}

static int siss_con_lslv(struct smtx *m, struct vec *x, struct vec *f)
{
    int n = m->pps.n;

    int *mia = m->ia;
    int *mja = m->ja;

    double *mdr = m->dr;
    double *mlr = m->lr;

    double *xv = x->dat;
    double *fv = f->dat;

    for (int i = 0; i < n; ++i) {
        int lr0 = mia[i];
        int lr1 = mia[i + 1];

        double s = 0;

        for (int lr = lr0; lr < lr1; ++lr) {
            int j = mja[lr];

            s += mlr[lr] * xv[j];
        }

        xv[i] = (fv[i] - s) / mdr[i];
    }

    return 0;
}

static int siss_con_uslv(struct smtx *m, struct vec *x, struct vec *f)
{
    int n = m->pps.n;

    int *mia = m->ia;
    int *mja = m->ja;

    double *mdr = m->dr;
    double *mur = m->ur;

    double *xv = x->dat;
    double *fv = f->dat;

    if (x != f)
        memcpy(xv, fv, sizeof(double) * n);

    for (int j = n - 1; j > -1; --j) {
        xv[j] = xv[j] / mdr[j];

        double x = xv[j];

        int ur0 = mia[j];
        int ur1 = mia[j + 1];

        for (int ur = ur0; ur < ur1; ++ur) {
            int ui = mja[ur];

            xv[ui] -= mur[ur] * x;
        }
    }

    return 0;
}

static int siss_bcg_con_slv(
    struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
{
    int n = m->pps.n;
    int c = 0;

    struct vec r;
    struct vec z;
    struct vec p;
    struct vec s;
    struct vec h;

    struct vec pt;
    struct vec st;

    if (vec_new(&r, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&z, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&p, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&s, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&h, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&pt, n)) {
        c = -1;
        goto end;
    }

    if (vec_new(&st, n)) {
        c = -1;
        goto end;
    }

    double alf = 0;
    double bet = 0;
    double omg = 0;

    double tmp = 0;
    double nrm = 0;

    mtx_vmul(m, x, &r);
    vec_cmb(f, &r, &r, -1);

    vec_dup(&r, &z);
    vec_dup(&r, &p);

    for (int k = 1; k <= o->ops.max; ++k) {
        siss_con_lslv(o->con.sm, &pt, &p);
        siss_con_uslv(o->con.sm, &pt, &pt);

        mtx_vmul(m, &pt, &h);

        vec_dot(&r, &z, &tmp);
        vec_dot(&h, &z, &bet);

        alf = tmp / bet;

        vec_cmb(&r, &h, &s, -alf);

        siss_con_lslv(o->con.sm, &st, &s);
        siss_con_uslv(o->con.sm, &st, &st);

        mtx_vmul(m, &s, &r);
        vec_dot(&r, &s, &omg);
        vec_dot(&r, &r, &bet);

        omg = omg / bet;

        vec_cmb(x, &pt, x, alf);
        vec_cmb(x, &st, x, omg);
        vec_cmb(&s, &r, &r, -omg);
        vec_nrm(&r, &nrm);

        if (o->ops.itr.run)
            o->ops.itr.run(o->ops.itr.ctx, &o->ops);

        if (nrm < o->ops.err)
            break;

        vec_dot(&r, &z, &bet);

        bet = (bet / tmp) * (alf / omg);

        vec_cmb(&p, &h, &p, -omg);
        vec_cmb(&r, &p, &p, bet);
    }

end:
    vec_cls(&r);
    vec_cls(&z);
    vec_cls(&p);
    vec_cls(&s);
    vec_cls(&h);

    vec_cls(&st);
    vec_cls(&pt);

    return c;
}

int siss_bcg_slv(
    struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
{
    if (!m || !x || !f) {
        errno = EINVAL;
        return -1;
    }

    return o->con.sm ? siss_bcg_con_slv(m, x, f, o)
                     : siss_bcg_unc_slv(m, x, f, o);
}

static void rot_inc(
    struct imtx *q, struct imtx *x, struct vec *b, double c, double s);

static void rot_mul_mtx(struct imtx *q, struct imtx *h, struct imtx *r);
static void rot_mul_vec(struct imtx *q, double b, struct vec *g);

static void red_slv(struct imtx *r, struct vec *y, struct vec *g);

int siss_gmr_slv(
    struct smtx *sm, struct vec *vx, struct vec *vf, struct iss_gmr_ops *ops)
{
    assert(sm);
    assert(vx);
    assert(vf);
    assert(ops);

    int m = ops->ops.max;
    int j = 0;

    double e = ops->ops.err;
    double n = vf->n;
    double b = 0;

    struct imtx v = {
        .pps = {.n = n, .m = 0},
        .dat = malloc(sizeof(double *) * m),
    };
    struct imtx h = {
        .pps = {.n = 1, .m = 0},
        .dat = malloc(sizeof(double *) * m),
    };
    struct imtx q = {
        .pps = {.n = 1, .m = 1},
        .dat = malloc(sizeof(double *) * (m + 1)),
    };
    struct imtx r = {
        .pps = {.n = 0, .m = 0},
        .dat = malloc(sizeof(double *) * m),
    };
    struct imtx x = {
        .pps = {.n = 1, .m = 1},
        .dat = malloc(sizeof(double *) * (m + 2)),
    };

    struct vec o;
    struct vec g;
    struct vec y;
    struct vec t;

    vec_new(&o, n);
    vec_new(&g, m);
    vec_new(&y, m);
    vec_new(&t, m + 1);

    mtx_vmul(sm, vx, &o);
    vec_cmb(vf, &o, &o, -1);
    vec_nrm(&o, &b);

    g.n = 1;
    y.n = 0;

    v.dat[0] = calloc(n, sizeof(double));
    q.dat[0] = calloc(m + 1, sizeof(double));
    x.dat[0] = calloc(m + 1, sizeof(double));

    for (int i = 0; i < n; ++i)
        v.dat[0][i] = o.dat[i] / b;

    for (j = 0; j < m; ++j) {
        h.dat[j] = calloc(j + 2, sizeof(double));
        r.dat[j] = calloc(j + 1, sizeof(double));
        q.dat[j + 1] = calloc(m + 1, sizeof(double));
        x.dat[j + 1] = calloc(m + 1, sizeof(double));

        v.pps.m += 1;
        h.pps.n += 1;
        h.pps.m += 1;
        q.pps.n += 1;
        q.pps.m += 1;
        x.pps.n += 1;
        x.pps.m += 1;
        r.pps.n += 1;
        r.pps.m += 1;

        g.n += 1;
        y.n += 1;

        struct vec vj = {.n = n, .dat = v.dat[j]};
        struct vec vi = {.n = n};
        struct vec hj = {.n = j + 2, .dat = h.dat[j]};

        mtx_vmul(sm, &vj, &o);

        for (int i = 0; i <= j; ++i) {
            vi.dat = v.dat[i];

            vec_dot(&o, &vi, &hj.dat[i]);
            vec_cmb(&o, &vi, &o, -hj.dat[i]);
        }

        vec_nrm(&o, &hj.dat[j + 1]);

        double hjj = hj.dat[j];
        double hjn = hj.dat[j + 1];
        double div = sqrt(hjj * hjj + hjn * hjn);

        double c = hjj / div;
        double s = hjn / div;

        rot_inc(&q, &x, &t, c, s);
        rot_mul_mtx(&q, &h, &r);
        rot_mul_vec(&q, b, &g);

        if (fabs(g.dat[j + 1]) < e)
            break;

        if (h.dat[j][j + 1] == 0 || j == m - 1)
            break;

        v.dat[j + 1] = calloc(n, sizeof(double));

        struct vec vn = {.n = n, .dat = v.dat[j + 1]};

        vec_mul(&o, &vn, 1.0 / hj.dat[j + 1]);
    }

    if (j == m)
        return 0;

    red_slv(&r, &y, &g);

    for (int i = 0; i < n; ++i)
        for (int p = 0; p < j; ++p)
            vx->dat[i] += v.dat[p][i] * y.dat[p];

    for (int i = 0; i < j; ++i) {
        free(v.dat[i]);
        free(h.dat[i]);
        free(q.dat[i]);
        free(r.dat[i]);
        free(x.dat[i]);
    }

    free(x.dat[j]);
    free(q.dat[j]);
    free(v.dat);
    free(h.dat);
    free(q.dat);
    free(r.dat);

    vec_cls(&o);
    vec_cls(&g);
    vec_cls(&y);
    vec_cls(&t);

    return 0;
}

static void rot_inc(
    struct imtx *q, struct imtx *x, struct vec *b, double c, double s)
{
    int n = q->pps.n;

    if (n == 2) {
        q->dat[0][0] = c;
        q->dat[0][1] = -s;
        q->dat[1][0] = s;
        q->dat[1][1] = c;

        return;
    }

    q->dat[n - 1][n - 1] = 1;

    x->dat[n - 2][n - 2] = c;
    x->dat[n - 1][n - 1] = c;
    x->dat[n - 1][n - 2] = -s;
    x->dat[n - 2][n - 1] = s;
    x->dat[n - 3][n - 2] = 0;
    x->dat[n - 2][n - 3] = 0;
    x->dat[n - 3][n - 3] = 1;

    struct vec ri = {.n = n};
    struct vec cj = {.n = n};

    for (int j = 0; j < n; ++j) {
        cj.dat = q->dat[j];

        for (int i = 0; i < n; ++i) {
            ri.dat = x->dat[i];
            vec_dot(&ri, &cj, &b->dat[i]);
        }

        vec_dup(b, &cj);
    }
}

static void rot_mul_mtx(struct imtx *q, struct imtx *h, struct imtx *r)
{
    int n = r->pps.n;

    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            double s = 0;

            for (int e = 0; e < n; ++e)
                s += q->dat[e][i] * h->dat[j][e];

            r->dat[j][i] = s;
        }
}

static void rot_mul_vec(struct imtx *q, double b, struct vec *g)
{
    struct vec v = {.n = g->n, .dat = q->dat[0]};
    vec_mul(&v, g, b);
}

static void red_slv(struct imtx *r, struct vec *y, struct vec *g)
{
    int n = r->pps.n;

    for (int i = n - 1; i >= 0; --i) {
        double s = g->dat[i];

        for (int j = i + 1; j < n; ++j)
            s -= y->dat[j] * r->dat[j][i];

        y->dat[i] = s / r->dat[i][i];
    }
}
