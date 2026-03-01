#include "fem.h"

static const int MU[8] = {0, 1, 0, 1, 0, 1, 0, 1};
static const int NU[8] = {0, 0, 1, 1, 0, 0, 1, 1};
static const int TT[8] = {0, 0, 0, 0, 1, 1, 1, 1};

static const double S[2][2] = {
    {1.0,  -1.0},
    {-1.0, 1.0 }
};

static const double M[2][2] = {
    {2.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 6.0}
};

static double sx[2][2];
static double sy[2][2];
static double sz[2][2];

static double mx[2][2];
static double my[2][2];
static double mz[2][2];

static int qud_dim(struct sim *sim, struct qud *qud, double *hxi, double *hzt);

int fem_hmc_lin_asm(struct sim *sim, struct fem_hmc_ctx *ctx)
{
    struct v3d *vtx = sim->msh->vtx.v3d.dat;
    double      frq = sim->ops.hmc.frq;

    struct sim_fun_ctx fun_ctx = {
        .sim = sim,
        .vtx = -1,
        .qud = -1,
        .hxd = -1,
    };

    for (int hi = 0; hi < sim->msh->hxd.len; ++hi) {
        struct hxd *hxd = &sim->msh->hxd.dat[hi];
        struct obj *obj = &sim->obj.dat[hxd->pid];
        struct mat *mat = &sim->mat.dat[obj->mat];
        struct val *src = &sim->src.dat[obj->src];

        double nlam = mat->lam.as.num;
        double nsig = mat->sig.as.num;
        double nchi = mat->chi.as.num;
        double ssrc[8];
        double csrc[8];

        for (int i = 0; i < 8; ++i) {
            fun_ctx.vtx = hxd->vtx[i];
            fun_ctx.hxd = hi;

            struct vec vw = {.dat = vtx[fun_ctx.vtx].dat, .n = 3};

            ssrc[i] = src->as.hmc.sin(&fun_ctx, &vw);
            csrc[i] = src->as.hmc.cos(&fun_ctx, &vw);
        }

        int v0 = hxd->vtx[0];
        int v7 = hxd->vtx[7];

        double hx = vtx[v7].dat[0] - vtx[v0].dat[0];
        double hy = vtx[v7].dat[1] - vtx[v0].dat[1];
        double hz = vtx[v7].dat[2] - vtx[v0].dat[2];

        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                sx[i][j] = S[i][j] / hx;
                sy[i][j] = S[i][j] / hy;
                sz[i][j] = S[i][j] / hz;

                mx[i][j] = M[i][j] * hx;
                my[i][j] = M[i][j] * hy;
                mz[i][j] = M[i][j] * hz;
            }

        for (int i = 0; i < 8; ++i) {
            int gli = hxd->vtx[i];
            int mui = MU[i];
            int nui = NU[i];
            int tti = TT[i];

            double bis = 0;
            double bic = 0;

            for (int j = 0; j < 8; ++j) {
                int glj = hxd->vtx[j];
                int muj = MU[j];
                int nuj = NU[j];
                int ttj = TT[j];

                double mij = mx[muj][mui] * my[nuj][nui] * mz[ttj][tti];
                double sij = sx[muj][mui] *
                             my[nuj][nui] *
                             mz[ttj][tti] +
                             mx[muj][mui] *
                             sy[nuj][nui] *
                             mz[ttj][tti] +
                             mx[muj][mui] *
                             my[nuj][nui] *
                             sz[ttj][tti];

                double pij = nlam * sij - frq * frq * nchi * mij;
                double cij = frq * nsig * mij;

                if (mtx_inc(&ctx->mtx, gli * 2, glj * 2, pij))
                    return -1;

                if (mtx_inc(&ctx->mtx, gli * 2 + 1, glj * 2 + 1, pij))
                    return -1;

                if (mtx_inc(&ctx->mtx, gli * 2, glj * 2 + 1, -cij))
                    return -1;

                if (mtx_inc(&ctx->mtx, gli * 2 + 1, glj * 2, cij))
                    return -1;

                bis += ssrc[j] * mij;
                bic += csrc[j] * mij;
            }

            ctx->vec.dat[gli * 2] += bis;
            ctx->vec.dat[gli * 2 + 1] += bic;
        }
    }

    for (int qi = 0; qi < sim->msh->qud.len; ++qi) {
        struct qud     *qud = &sim->msh->qud.dat[qi];
        struct bnd     *bnd = &sim->bnd.dat[qud->pid];
        struct cnd_bnd *cnd = &sim->cnd_bnd.dat[bnd->cnd];

        double hxi = 0;
        double hzt = 0;

        qud_dim(sim, qud, &hxi, &hzt);

        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                mx[i][j] = M[i][j] * hxi;
                mz[i][j] = M[i][j] * hzt;
            }

        switch (cnd->type) {
            case CND_BND_NEU: {
                struct val *tta = &cnd->pps.neu.tta;

                double stta[4];
                double ctta[4];

                for (int i = 0; i < 4; ++i) {
                    fun_ctx.vtx = qud->vtx[i];
                    fun_ctx.qud = qi;
                    fun_ctx.hxd = qud->hxd;

                    struct vec vw = {.dat = vtx[fun_ctx.vtx].dat, .n = 3};

                    stta[i] = tta->as.hmc.sin(&fun_ctx, &vw);
                    ctta[i] = tta->as.hmc.cos(&fun_ctx, &vw);
                }

                for (int i = 0; i < 4; ++i) {
                    int gli = qud->vtx[i];
                    int mui = MU[i];
                    int nui = NU[i];

                    double bis = 0;
                    double bic = 0;

                    for (int j = 0; j < 4; ++j) {
                        int muj = MU[j];
                        int nuj = NU[j];

                        double mij = mx[muj][mui] * mz[nuj][nui];

                        bis += stta[j] * mij;
                        bic += ctta[j] * mij;
                    }

                    ctx->vec.dat[gli * 2] += bis;
                    ctx->vec.dat[gli * 2 + 1] += bic;
                }

                break;
            }
            case CND_BND_ROB: {
                struct val *bet = &cnd->pps.rob.bet;
                struct val *ext = &cnd->pps.rob.ext;

                double nbet = bet->as.num;
                double sext[4];
                double cext[4];

                for (int i = 0; i < 4; ++i) {
                    fun_ctx.vtx = qud->vtx[i];
                    fun_ctx.qud = qi;
                    fun_ctx.hxd = qud->hxd;

                    struct vec vw = {.dat = vtx[fun_ctx.vtx].dat, .n = 3};

                    sext[i] = ext->as.hmc.sin(&fun_ctx, &vw);
                    cext[i] = ext->as.hmc.cos(&fun_ctx, &vw);
                }

                for (int i = 0; i < 4; ++i) {
                    int gli = qud->vtx[i];
                    int mui = MU[i];
                    int nui = NU[i];

                    double bis = 0;
                    double bic = 0;

                    for (int j = 0; j < 4; ++j) {
                        int glj = qud->vtx[j];
                        int muj = MU[j];
                        int nuj = NU[j];

                        double mij = mx[muj][mui] * mz[nuj][nui];
                        double pij = nbet * mij;

                        if (mtx_inc(&ctx->mtx, gli * 2, glj * 2, pij))
                            return -1;

                        if (mtx_inc(&ctx->mtx, gli * 2 + 1, glj * 2 + 1, pij))
                            return -1;

                        bis += sext[j] * mij;
                        bic += cext[j] * mij;
                    }

                    ctx->vec.dat[gli * 2] += bis * nbet;
                    ctx->vec.dat[gli * 2 + 1] += bic * nbet;
                }

                break;
            }
            default:
                errno = ENOTSUP;
                return -1;
        }
    }

    return 0;
}

static int qud_dim(struct sim *sim, struct qud *qud, double *hxi, double *hzt)
{
    double dat[3];

    struct vec nrm = {
        .n = 3,
        .dat = dat,
    };

    if (umsh_qud_nrm(sim->msh, qud, &nrm))
        return -1;

    struct v3d *vtx = sim->msh->vtx.v3d.dat;

    int v0 = qud->vtx[0];
    int v3 = qud->vtx[3];

    if (nrm.dat[0] != 0) {
        *hxi = vtx[v3].dat[1] - vtx[v0].dat[1];
        *hzt = vtx[v3].dat[2] - vtx[v0].dat[2];
    } else if (nrm.dat[1] != 0) {
        *hxi = vtx[v3].dat[0] - vtx[v0].dat[0];
        *hzt = vtx[v3].dat[2] - vtx[v0].dat[2];
    } else {
        *hxi = vtx[v3].dat[0] - vtx[v0].dat[0];
        *hzt = vtx[v3].dat[1] - vtx[v0].dat[1];
    }

    return 0;
}
