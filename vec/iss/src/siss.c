#include <errno.h>
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

int siss_gmr_slv(
    struct smtx *m, struct vec *x, struct vec *f, struct iss_gmr_ops *ops)
{
    assert(m);
    assert(x);
    assert(f);
    assert(ops);

    double max = ops->ops.max;
    double n = f->n;
    double b = 0;

    struct imtx v;
    struct imtx h;
    struct imtx q;

    v.pps.n = n;
    v.pps.m = 1;
    v.dat = malloc(sizeof(double *) * max);
    v.dat[0] = malloc(sizeof(double) * n);

    h.pps.n = 1;
    h.pps.m = 0;
    h.dat = malloc(sizeof(double *) * max);

    q.pps.n = 0;
    q.pps.m = 0;
    q.dat = malloc(sizeof(double *) * max);

    memset(v.dat, 0, sizeof(double *) * n);
    memset(h.dat, 0, sizeof(double *) * n);

    struct vec r0;
    struct vec om;

    vec_new(&r0, n);
    vec_new(&om, n);

    mtx_vmul(m, x, &r0);
    vec_cmb(f, &r0, &r0, -1);
    vec_nrm(&r0, &b);

    for (int i = 0; i < n; ++i) {
        v.dat[0][i] = r0.dat[i] / b;
        q.dat = malloc(sizeof(double) * n);

        memset(q.dat, 0, sizeof(double) * n);
    }

    for (int j = 0; j < max; ++j) {
        v.pps.m += 1;
        h.pps.n += 1;
        h.pps.m += 1;

        struct vec vj = {.n = n, .dat = v.dat[j]};
        struct vec vi = {.n = n};
        struct vec vn;
        struct vec hj;

        vec_new(&vn, n);
        vec_new(&hj, j + 2);

        v.dat[j + 1] = vn.dat;
        h.dat[j] = hj.dat;

        mtx_vmul(m, &vj, &om);

        for (int i = 0; i <= j; ++i) {
            vi.dat = v.dat[i];

            vec_dot(&om, &vi, &hj.dat[i]);
            vec_cmb(&om, &vi, &om, -hj.dat[i]);
        }

        vec_nrm(&om, &hj.dat[j + 1]);

        if (h.dat[j][j + 1] == 0)
            break;

        vec_mul(&om, &vn, 1.0 / hj.dat[j + 1]);
    }

    return 0;
}
