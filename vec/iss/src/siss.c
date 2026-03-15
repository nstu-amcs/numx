#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <numx/vec/iss.h>

static int siss_bcg_unc_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
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

        if (o->ops.itr.run) {
            o->ops.itr.run(o->ops.itr.ctx, &o->ops);
        }

        if (nrm < o->ops.err) {
            break;
        }

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

static int siss_bcg_con_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
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

    vec_cls(&st);
    vec_cls(&pt);

    return c;
}

int siss_bcg_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *o)
{
    if (!m || !x || !f) {
        errno = EINVAL;
        return -1;
    }

    int r = o->con.sm ? siss_bcg_con_slv(m, x, f, o) : siss_bcg_unc_slv(m, x, f, o);

    printf("[vec][iss][bcg] itr: %d, err: %.7e\n", o->ops.run.itr, o->ops.run.err);

    return r;
}

static void red_slv(struct imtx *r, struct vec *y, struct vec *g);

int siss_gmr_slv(struct smtx *sm, struct vec *vx, struct vec *vf, struct iss_gmr_ops *ops)
{
    assert(sm);
    assert(vx);
    assert(vf);
    assert(ops);

    int m = ops->ops.max;
    int n = vf->n;
    int j = 0;

    double e = ops->ops.err;
    double b = 0;

    struct imtx v = {
        .pps = {.c = n, .r = 1},
        .dat = malloc(sizeof(double *) * m),
    };

    struct imtx h = {
        .pps = {.c = 1, .r = 0},
        .dat = malloc(sizeof(double *) * m),
    };

    struct vec o; // omega
    struct vec g; // right-hand side
    struct vec t; // buffer
    struct vec c; // rotation coefficients

    vec_new(&o, n);
    vec_new(&g, m + 1);
    vec_new(&t, m + 1);
    vec_new(&c, m * 2);

    mtx_vmul(sm, vx, &o);
    vec_cmb(vf, &o, &o, -1);
    vec_nrm(&o, &b);

    if (b < e) {
        ops->ops.run.itr = 0;
        ops->ops.run.err = b;

        if (ops->ops.itr.run) {
            ops->ops.itr.run(ops->ops.itr.ctx, &ops->ops);
        }

        goto end;
    }

    g.n = 1;
    v.dat[0] = malloc(sizeof(double) * n);
    g.dat[0] = b;

    for (int i = 0; i < n; ++i) {
        v.dat[0][i] = o.dat[i] / b;
    }

    for (j = 0; j < m; ++j) {
        h.dat[j] = malloc(sizeof(double) * (j + 2));
        h.pps.c += 1;
        h.pps.r += 1;
        g.n += 1;

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

        double hsv = h.dat[j][j + 1];

        for (int i = 0; i < j; ++i) {
            hj.n = i + 2;
            t.n = i + 2;

            vec_rot(&hj, &t, i, c.dat[i * 2], c.dat[i * 2 + 1]);
            vec_dup(&t, &hj);
        }

        double hjj = hj.dat[j];
        double hjn = hj.dat[j + 1];
        double div = sqrt(hjj * hjj + hjn * hjn);

        if (div == 0.0) {
            break;
        }

        c.dat[j * 2] = hjj / div;
        c.dat[j * 2 + 1] = hjn / div;

        hj.n = j + 2;
        t.n = j + 2;

        vec_rot(&hj, &t, j, c.dat[j * 2], c.dat[j * 2 + 1]);
        vec_dup(&t, &hj);

        vec_rot(&g, &t, j, c.dat[j * 2], c.dat[j * 2 + 1]);
        vec_dup(&t, &g);

        ops->ops.run.itr = j;
        ops->ops.run.err = fabs(g.dat[j + 1]);

        if (ops->ops.itr.run) {
            ops->ops.itr.run(ops->ops.itr.ctx, &ops->ops);
        }

        if (hsv == 0 || fabs(g.dat[j + 1]) < e) {
            break;
        }

        if (j < m - 1) {
            v.dat[j + 1] = malloc(sizeof(double) * n);
            v.pps.r += 1;

            struct vec vn = {.n = n, .dat = v.dat[j + 1]};
            vec_mul(&o, &vn, 1.0 / hsv);
        }
    }

    if (j == m) {
        j -= 1;
    }

    h.pps.c = j + 1;
    o.n = j + 1;
    g.n = j + 1;

    red_slv(&h, &o, &g);

    for (int i = 0; i < n; ++i)
        for (int p = 0; p <= j; ++p)
            vx->dat[i] += v.dat[p][i] * o.dat[p];

end:
    printf("[vec][iss][gmr] itr: %d, err: %.7e\n", ops->ops.run.itr, ops->ops.run.err);

    for (int i = 0; i < h.pps.r; ++i) {
        free(h.dat[i]);
    }

    for (int i = 0; i < v.pps.r; ++i) {
        free(v.dat[i]);
    }

    free(v.dat);
    free(h.dat);

    vec_cls(&o);
    vec_cls(&g);
    vec_cls(&t);
    vec_cls(&c);

    return 0;
}

static void red_slv(struct imtx *r, struct vec *y, struct vec *g)
{
    int n = r->pps.c;

    for (int i = n - 1; i >= 0; --i) {
        double s = g->dat[i];

        for (int j = i + 1; j < n; ++j)
            s -= y->dat[j] * r->dat[j][i];

        y->dat[i] = s / r->dat[i][i];
    }
}
