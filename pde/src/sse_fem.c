#include <assert.h>
#include <numx/pde/sse.h>
#include <numx/vec/iss.h>
#include <numx/vec/mtx.h>
#include <stdbool.h>

const double G[2][2] = {
    {1.0,  -1.0},
    {-1.0, 1.0 }
};

const double M[2][2] = {
    {2.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 6.0}
};

const double X[2][2] = {
    {-0.5, 0.5},
    {-0.5, 0.5}
};

const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

const int DEM[8][8] = {
    {0,  1,  3,  4,  9,  10, 12, 13},
    {14, 0,  2,  3,  8,  9,  11, 12},
    {16, 15, 0,  1,  6,  7,  9,  10},
    {17, 16, 14, 0,  5,  6,  8,  9 },
    {22, 21, 19, 18, 0,  1,  3,  4 },
    {23, 22, 20, 19, 14, 0,  2,  3 },
    {25, 24, 22, 21, 16, 15, 0,  1 },
    {26, 25, 23, 22, 17, 16, 14, 0 }
};

const int DFM[3][4][4] = {
    {{0, 1, 3, 4},  {14, 0, 2, 3}, {16, 15, 0, 1}, {17, 16, 14, 0}},
    {{0, 3, 9, 12}, {16, 0, 6, 9}, {22, 19, 0, 3}, {26, 22, 16, 0}},
    {{0, 1, 9, 10}, {14, 0, 8, 9}, {22, 21, 0, 1}, {23, 22, 14, 0}}
};

static int hxd_evo(struct obj *obj, struct dmtx *m, struct vec *f, bool *stat);

static int qud_evo_neu(struct obj *obj, struct qud *qud, struct vec *f);
static int qud_evo_rob(struct obj *obj, struct qud *qud, struct dmtx *m, struct vec *f);
static int qud_evo_dir(struct obj *obj, struct qud *qud, struct dmtx *m, struct vec *f);

int pde_sse_fem_slv(struct obj *obj, struct vec *x, struct sse_fem_ops ops)
{
    assert(obj);
    assert(x);

    int r = 0;
    int n = obj->ax.len * obj->ay.len * obj->az.len;

    struct dmtx m;
    struct vec  f;

    bool *stat = malloc(sizeof(bool) * n);

    if (!stat)
        goto end;

    memset(stat, 0, sizeof(bool) * n);

    if (mtx_new(&m, ((struct dmtx_pps){.n = n, .d = 27})))
        goto end;

    for (int i = 0; i < n; ++i)
        m.ad[i][0] = 1;

    if (vec_new(&f, n))
        goto end;

    if (hxd_evo(obj, &m, &f, stat))
        goto end;

    struct qud *qud = obj->qud.dat;

    for (int q = 0; q < obj->qud.len; ++q) {
        struct fctx *ctx = qud[q].ctx;

        switch (ctx->cnd.mod) {
            case CND_DIR:
                continue;
            case CND_NEU:
                if (qud_evo_neu(obj, &qud[q], &f))
                    goto end;

                break;
            case CND_ROB:
                if (qud_evo_rob(obj, &qud[q], &m, &f))
                    goto end;

                break;
        }
    }

    for (int q = 0; q < obj->qud.len; ++q) {
        struct fctx *ctx = qud[q].ctx;

        if (ctx->cnd.mod == CND_DIR)
            if (qud_evo_dir(obj, &qud[q], &m, &f))
                goto end;
    }

    for (int i = 0; i < n; ++i)
        if (!stat[i])
            m.ad[i][0] = 1;

    if (vec_new(x, n))
        goto end;

    switch (ops.ops.iss.mod) {
        case ISS_JAC:
            if ((r = iss_jac_slv(&m, x, &f, ops.ops.iss.ops.jac)))
                goto end;

            break;
        case ISS_RLX:
            if ((r = iss_rlx_slv(&m, x, &f, ops.ops.iss.ops.rlx)))
                goto end;

            break;

        default:
            r = -1;
            errno = ENOTSUP;
            goto end;
    }

end:
    mtx_cls(&m);
    vec_cls(&f);

    free(stat);

    return 0;
}

static int hxd_evo(struct obj *obj, struct dmtx *m, struct vec *f, bool *stat)
{
    assert(obj);
    assert(m);
    assert(f);

    int nx = obj->ax.len;
    int ny = obj->ay.len;

    double *ax = obj->ax.dat;
    double *ay = obj->ay.dat;
    double *az = obj->az.dat;

    int map[8] = {0, 1, nx, nx + 1, nx * ny, nx * ny + 1, nx * (ny + 1), nx * (ny + 1) + 1};

    for (int h = 0; h < obj->hxd.len; ++h) {
        struct ectx *ctx = obj->hxd.dat[h].ctx;

        double lam = ctx->lam;
        double gam = ctx->gam;

        int v0 = obj->hxd.dat[h].vtx[0];
        int v1 = obj->hxd.dat[h].vtx[1];
        int v2 = obj->hxd.dat[h].vtx[2];
        int v4 = obj->hxd.dat[h].vtx[4];

        int x0 = v0 % nx;
        int x1 = v1 % nx;
        int y0 = v0 / nx;
        int y1 = v2 / nx;
        int z0 = v0 / (nx * ny);
        int z1 = v4 / (nx * ny);

        double gx[2][2];
        double gy[2][2];
        double gz[2][2];

        double mx[2][2];
        double my[2][2];
        double mz[2][2];

        int dx = nx - (x1 - x0);
        int dy = (nx * ny) - ((x1 - x0) * (y1 - y0));

        int n = v0;

        for (int zi = z0; zi < z1; ++zi, n += dy)
            for (int yi = y0; yi < y1; ++yi, n += dx)
                for (int xi = x0; xi < x1; ++xi, n += 1) {
                    double vx[2] = {ax[xi], ax[xi + 1]};
                    double vy[2] = {ay[yi], ay[yi + 1]};
                    double vz[2] = {az[zi], az[zi + 1]};

                    double hx = vx[1] - vx[0];
                    double hy = vy[1] - vy[0];
                    double hz = vz[1] - vz[0];

                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j) {
                            gx[i][j] = G[i][j] / hx;
                            gy[i][j] = G[i][j] / hy;
                            gz[i][j] = G[i][j] / hz;

                            mx[i][j] = M[i][j] * hx;
                            my[i][j] = M[i][j] * hy;
                            mz[i][j] = M[i][j] * hz;
                        }

                    double dec[8];

                    for (int i = 0; i < 8; ++i)
                        dec[i] = ctx->ext(&((struct vtx){.x = vx[MU[i]], .y = vy[NU[i]], .z = vz[TT[i]]}));

                    for (int i = 0; i < 8; ++i) {
                        int gi = n + map[i];

                        stat[gi] = true;

                        int mui = MU[i];
                        int nui = NU[i];
                        int tti = TT[i];

                        for (int j = 0; j < 8; ++j) {
                            int muj = MU[j];
                            int nuj = NU[j];
                            int ttj = TT[j];

                            double mxl = mx[mui][muj];
                            double myl = my[nui][nuj];
                            double mzl = mz[tti][ttj];

                            double gxl = gx[mui][muj];
                            double gyl = gy[nui][nuj];
                            double gzl = gz[tti][ttj];

                            m->ad[gi][DEM[i][j]] += gam * mxl * myl * mzl;
                            m->ad[gi][DEM[i][j]] += lam * (gxl * myl * mzl + mxl * gyl * mzl + mxl * myl * gzl);

                            f->dat[gi] += dec[j] * mxl * myl * mzl;
                        }
                    }
                }
    }

    return 0;
}

static int qud_evo_neu(struct obj *obj, struct qud *qud, struct vec *f)
{
    struct fctx *ctx = qud->ctx;

    int nx = obj->ax.len;
    int ny = obj->ay.len;

    double *ax = obj->ax.dat;
    double *ay = obj->ay.dat;
    double *az = obj->az.dat;

    int v0 = qud->vtx[0];
    int v1 = qud->vtx[1];
    int v2 = qud->vtx[2];

    double *axi = 0;
    double *azt = 0;

    struct vtx v = {.ctx = 0, .x = 0, .y = 0, .z = 0, .n = 0};

    double *xi = 0;
    double *zt = 0;

    int xi0 = 0;
    int xi1 = 0;
    int xis = 0;

    int zt0 = 0;
    int zt1 = 0;
    int zts = 0;

    int map[4];

    enum norm nrm = qud_norm(obj, qud);

    switch (nrm) {
        case NORM_L:
        case NORM_R:
            axi = ay;
            azt = az;

            v.x = ax[v0 % nx];
            xi = &v.y;
            zt = &v.z;

            xi0 = v0 / nx;
            xi1 = v1 / nx;
            xis = nx;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            map[0] = 0;
            map[1] = nx;
            map[2] = nx * ny;
            map[3] = nx * ny + nx;

            break;
        case NORM_F:
        case NORM_B:
            axi = ax;
            azt = az;

            xi = &v.x;
            v.y = ay[v0 / nx];
            zt = &v.z;

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            map[0] = 0;
            map[1] = 1;
            map[2] = nx * ny;
            map[3] = nx * ny + 1;

            break;
        case NORM_U:
        case NORM_D:
            axi = ax;
            azt = ay;

            xi = &v.x;
            zt = &v.y;
            v.z = az[v0 / (nx * ny)];

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / nx;
            zt1 = v2 / nx;
            zts = nx;

            map[0] = 0;
            map[1] = 1;
            map[2] = nx;
            map[3] = nx + 1;

            break;
    }

    double hxi = 0;
    double hzt = 0;

    double mxi[2][2];
    double mzt[2][2];
    double dec[4];

    int n = v0;

    for (int xii = xi0; xii < xi1; ++xii) {
        int nxi = n;

        hxi = axi[xii + 1] - axi[xii];

        for (int zti = zt0; zti < zt1; ++zti) {
            hzt = azt[zti + 1] - azt[zti];

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) {
                    mxi[i][j] = hxi * M[i][j];
                    mzt[i][j] = hzt * M[i][j];
                }

            *xi = axi[xii];
            *zt = azt[zti];
            dec[0] = ctx->cnd.pps.neu.tta(&v);

            *xi = axi[xii + 1];
            *zt = azt[zti];
            dec[1] = ctx->cnd.pps.neu.tta(&v);

            *xi = axi[xii];
            *zt = azt[zti + 1];
            dec[2] = ctx->cnd.pps.neu.tta(&v);

            *xi = axi[xii + 1];
            *zt = azt[zti + 1];
            dec[3] = ctx->cnd.pps.neu.tta(&v);

            for (int i = 0; i < 4; ++i) {
                int mui = MU[i];
                int nui = NU[i];

                for (int j = 0; j < 4; ++j) {
                    int muj = MU[j];
                    int nuj = NU[j];

                    f->dat[n + map[i]] += dec[j] * mxi[mui][muj] * mzt[nui][nuj];
                }
            }

            n += zts;
        }

        n = nxi + xis;
    }

    return 0;
}

static int qud_evo_rob(struct obj *obj, struct qud *qud, struct dmtx *m, struct vec *f)
{
    struct fctx *ctx = qud->ctx;

    int nx = obj->ax.len;
    int ny = obj->ay.len;

    double *ax = obj->ax.dat;
    double *ay = obj->ay.dat;
    double *az = obj->az.dat;

    int v0 = qud->vtx[0];
    int v1 = qud->vtx[1];
    int v2 = qud->vtx[2];

    double *axi = 0;
    double *azt = 0;

    struct vtx v = {.ctx = 0, .x = 0, .y = 0, .z = 0, .n = 0};

    double *xi = 0;
    double *zt = 0;

    int xi0 = 0;
    int xi1 = 0;
    int xis = 0;

    int zt0 = 0;
    int zt1 = 0;
    int zts = 0;

    int map[4];

    enum norm nrm = qud_norm(obj, qud);

    switch (nrm) {
        case NORM_L:
        case NORM_R:
            axi = ay;
            azt = az;

            v.x = ax[v0 % nx];
            xi = &v.y;
            zt = &v.z;

            xi0 = v0 / nx;
            xi1 = v1 / nx;
            xis = nx;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            map[0] = 0;
            map[1] = nx;
            map[2] = nx * ny;
            map[3] = nx * ny + nx;

            break;
        case NORM_F:
        case NORM_B:
            axi = ax;
            azt = az;

            xi = &v.x;
            v.y = ay[v0 / nx];
            zt = &v.z;

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            map[0] = 0;
            map[1] = 1;
            map[2] = nx * ny;
            map[3] = nx * ny + 1;

            break;
        case NORM_U:
        case NORM_D:
            axi = ax;
            azt = ay;

            xi = &v.x;
            zt = &v.y;
            v.z = az[v0 / (nx * ny)];

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / nx;
            zt1 = v2 / nx;
            zts = nx;

            map[0] = 0;
            map[1] = 1;
            map[2] = nx;
            map[3] = nx + 1;

            break;
    }

    double hxi = 0;
    double hzt = 0;

    double mxi[2][2];
    double mzt[2][2];
    double dec[4];

    int n = v0;

    for (int xii = xi0; xii < xi1; ++xii) {
        int nxi = n;

        hxi = axi[xii + 1] - axi[xii];

        for (int zti = zt0; zti < zt1; ++zti) {
            hzt = azt[zti + 1] - azt[zti];

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) {
                    mxi[i][j] = hxi * M[i][j];
                    mzt[i][j] = hzt * M[i][j];
                }

            *xi = axi[xii];
            *zt = azt[zti];
            dec[0] = ctx->cnd.pps.rob.tmp(&v);

            *xi = axi[xii + 1];
            *zt = azt[zti];
            dec[1] = ctx->cnd.pps.rob.tmp(&v);

            *xi = axi[xii];
            *zt = azt[zti + 1];
            dec[2] = ctx->cnd.pps.rob.tmp(&v);

            *xi = axi[xii + 1];
            *zt = azt[zti + 1];
            dec[3] = ctx->cnd.pps.rob.tmp(&v);

            for (int i = 0; i < 4; ++i) {
                int mui = MU[i];
                int nui = NU[i];

                double fi = 0;

                for (int j = 0; j < 4; ++j) {
                    int muj = MU[j];
                    int nuj = NU[j];

                    m->ad[n + map[i]][DFM[nrm / 2][i][j]] += ctx->cnd.pps.rob.bet * mxi[mui][muj] * mzt[nui][nuj];
                    fi += dec[j] * mxi[mui][muj] * mzt[nui][nuj];
                }

                f->dat[n + map[i]] += fi * ctx->cnd.pps.rob.bet;
            }

            n += zts;
        }

        n = nxi + xis;
    }

    return 0;
}

static int qud_evo_dir(struct obj *obj, struct qud *qud, struct dmtx *m, struct vec *f)
{
    struct fctx *ctx = qud->ctx;

    int nx = obj->ax.len;
    int ny = obj->ay.len;

    double *ax = obj->ax.dat;
    double *ay = obj->ay.dat;
    double *az = obj->az.dat;

    int v0 = qud->vtx[0];
    int v1 = qud->vtx[1];
    int v2 = qud->vtx[2];

    double *axi = 0;
    double *azt = 0;

    struct vtx v = {.ctx = 0, .x = 0, .y = 0, .z = 0, .n = 0};

    double *xi = 0;
    double *zt = 0;

    int xi0 = 0;
    int xi1 = 0;
    int xis = 0;

    int zt0 = 0;
    int zt1 = 0;
    int zts = 0;

    enum norm nrm = qud_norm(obj, qud);

    switch (nrm) {
        case NORM_L:
        case NORM_R:
            axi = ay;
            azt = az;

            v.x = ax[v0 % nx];
            xi = &v.y;
            zt = &v.z;

            xi0 = v0 / nx;
            xi1 = v1 / nx;
            xis = nx;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            break;
        case NORM_F:
        case NORM_B:
            axi = ax;
            azt = az;

            xi = &v.x;
            v.y = ay[v0 / nx];
            zt = &v.z;

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / (nx * ny);
            zt1 = v2 / (nx * ny);
            zts = nx * ny;

            break;
        case NORM_U:
        case NORM_D:
            axi = ax;
            azt = ay;

            xi = &v.x;
            zt = &v.y;
            v.z = az[v0 / (nx * ny)];

            xi0 = v0 % nx;
            xi1 = v1 % nx;
            xis = 1;

            zt0 = v0 / nx;
            zt1 = v2 / nx;
            zts = nx;

            break;
    }

    int n = v0;

    for (int xii = xi0; xii <= xi1; ++xii) {
        int nxi = n;

        for (int zti = zt0; zti <= zt1; ++zti) {
            memset(m->ad[n], 0, sizeof(double) * 27);

            *xi = axi[xii];
            *zt = azt[zti];

            f->dat[n] = ctx->cnd.pps.dir.tmp(&v);
            m->ad[n][0] = 1;

            n += zts;
        }

        n = nxi + xis;
    }

    return 0;
}
