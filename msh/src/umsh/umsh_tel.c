#include <assert.h>
#include <stdarg.h>
#include <stdio.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/msh/umsh.h>

#define less(a, b) ((b) - (a) > 1e-7)

typedef struct area
{
    int pid;

    double x0;
    double x1;
    double y0;
    double y1;
} area;

log_def(area_log, area_rec, area);
log_gen(area_log, area_rec, area, PUB);

struct div_ops
{
    double ini;
    double inc;
    int    dir;
};

static int axis_fin(struct dlog *axis, struct div_ops *ops);
static int axis_div(struct dlog *axis, int f);

int umsh_imp_tel(struct umsh *msh, const char *dir, const char *pfx)
{
    assert(msh);
    assert(dir);
    assert(pfx);
    assert(msh->type == MSH_C2D);

    double t = 0;
    int    r = 0;
    int    n = 0;
    char   p[256];
    FILE  *f;

    snprintf(p, sizeof(p), "%s/%s.tel", dir, pfx);

    if (!(f = fopen(p, "r"))) {
        r = -1;
        goto end;
    }

    if (fscanf(f, "%d", &n) != 1) {
        r = -1;
        goto end;
    }

    area_log areas;

    if ((r = area_log_new(&areas)) == -1) {
        goto end;
    }

    struct area a;

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf %lf %lf %lf %lf %lf %d", &a.x0, &a.x1, &a.y0, &a.y1, &t, &t, &a.pid);
        a.pid -= 1;
        area_log_add(&areas, a);
    }

    dlog x;
    dlog y;

    if ((r = dlog_new(&x)) == -1) {
        goto end;
    }

    if ((r = dlog_new(&y)) == -1) {
        goto end;
    }

    fscanf(f, "%lf %d", &t, &n);
    dlog_add(&x, t);

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &t);
        dlog_add(&x, t);
    }

    struct div_ops *x_ops = malloc(sizeof(struct div_ops) * n);

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &x_ops[i].ini);
    }

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &x_ops[i].inc);
    }

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%d", &x_ops[i].dir);
    }

    fscanf(f, "%lf %d", &t, &n);
    dlog_add(&y, t);

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &t);
        dlog_add(&y, t);
    }

    struct div_ops *y_ops = malloc(sizeof(struct div_ops) * n);

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &y_ops[i].ini);
    }

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &y_ops[i].inc);
    }

    for (int i = 0; i < n; ++i) {
        fscanf(f, "%d", &y_ops[i].dir);
    }

    axis_fin(&x, x_ops);
    axis_fin(&y, y_ops);

    int sx, sy;

    fscanf(f, "%d %d", &sx, &sy);
    axis_div(&x, sx);
    axis_div(&y, sy);

    v2d_cut_dev(&msh->vtx.v2d, x.len * y.len);
    qud_cut_exp(&msh->qud, (x.len - 1) * (y.len - 1));
    dlog_rst(&y);

    double xv = 0;
    double yv = 0;
    int    i = 0;

    while (dlog_adv(&y, &yv) != -1) {
        dlog_rst(&x);

        while (dlog_adv(&x, &xv) != -1) {
            struct v2d *vtx = &msh->vtx.v2d.dat[i++];

            vtx->dat[0] = xv;
            vtx->dat[1] = yv;
        }
    }

    area_log_rst(&areas);
    int rs = x.len;
    int rp = 0;

    while (area_log_adv(&areas, &a) != -1) {
        dlog_rst(&y);
        dlog_adv(&y, &yv);
        rp = 0;

        while (less(yv, a.y0)) {
            dlog_adv(&y, &yv);
            rp += rs;
        }

        while (less(yv, a.y1)) {
            dlog_rst(&x);
            dlog_adv(&x, &xv);
            i = rp;

            while (less(xv, a.x0)) {
                dlog_adv(&x, &xv);
                i += 1;
            }

            while (less(xv, a.x1)) {
                struct qud qud = {
                    .pid = a.pid,
                };

                qud.vtx[0] = i;
                qud.vtx[1] = i + 1;
                qud.vtx[2] = i + rs;
                qud.vtx[3] = i + rs + 1;
                qud_cut_add(&msh->qud, qud);

                if (dlog_adv(&x, &xv) == -1) {
                    break;
                }

                i += 1;
            }

            if (dlog_adv(&y, &yv) == -1) {
                break;
            }

            rp += rs;
        }
    }

    printf("Areas: %d\n", areas.len);
    printf("Vertices: %d\n", msh->vtx.v2d.len);
    printf("Elements: %d\n", msh->qud.len);

end:
    dlog_cls(&x);
    dlog_cls(&y);
    area_log_cls(&areas);
    fclose(f);

    return r;
}

static int axis_fin(struct dlog *axis, struct div_ops *ops)
{
    double a = 0;
    double b = 0;
    double x = 0;
    double s = 0;
    int    i = 0;

    dlog_rst(axis);
    dlog_adv(axis, &a);

    while (dlog_adv(axis, &b) != -1) {
        s = ops[i].ini;
        x = a + s;

        while (less(x, b)) {
            dlog_ins(axis, x, LOG_L);
            x += s;
        }

        a = b;
        i += 1;
    }

    return 0;
}

static int axis_div(struct dlog *axis, int f)
{
    (void)axis;
    (void)f;

    double a = 0;
    double b = 0;
    double s = 0;
    double x = 0;

    dlog_rst(axis);
    dlog_adv(axis, &a);

    while (dlog_adv(axis, &b) != -1) {
        s = (b - a) / f;
        x = a + s;

        while (less(x, b)) {
            dlog_ins(axis, x, LOG_L);
            x += s;
        }

        a = b;
    }

    return 0;
}