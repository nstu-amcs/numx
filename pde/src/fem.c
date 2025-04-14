#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>

#include "fem_lin.h"

int fem_new(struct fem *fem)
{
    assert(fem);

    fem->slv.exe.run = fem_exe;

    fem->ops.mod = FEM_STD;
    fem->ops.bss = FEM_BSS_LIN;
    fem->ops.iss.mod = ISS_BCG;
    fem->ops.iss.ops.bcg.con.sm = 0;
    fem->ops.iss.ops.bcg.ops.err = 1e-10;
    fem->ops.iss.ops.bcg.ops.itr.run = 0;
    fem->ops.iss.ops.bcg.ops.max = 500;
    fem->ops.non.err = 1e-10;
    fem->ops.non.max = 50;
    fem->ops.non.rlx = 1;
    fem->ops.tdd = FEM_TDD_I2S;

    return 0;
}

int fem_cls(struct fem *fem)
{
    assert(fem);

    mtx_cls(&fem->prv.ell);
    mtx_cls(&fem->prv.pbc);
    mtx_cls(&fem->prv.hyp);
    vec_cls(&fem->prv.vec);

    return 0;
}

static int fem_ini(struct sim *sim);

int fem_exe(void *ctx, struct sim *sim)
{
    (void)ctx;

    if (fem_ini(sim))
        return -1;

    sim->slv->apx.ctx = 0;
    sim->slv->apx.run = fem_lin_apx;

    switch (sim->mod) {
        case SIM_ELL:
            return fem_lin_ell_slv(sim);
        case SIM_PBC:
            return fem_lin_pbc_slv(sim);
        case SIM_HYP:
            return fem_lin_hyp_slv(sim);
    }

    return 0;
}

static int fem_ini(struct sim *sim)
{
    struct fem *fem = (struct fem *)sim->slv;

    int n = sim->msh->vtx.len;
    int z = 0;
    int r = 0;

    struct ilog *map = malloc(sizeof(struct ilog) * n);

    if (!map) {
        r = -1;
        goto end;
    }

    for (int i = 0; i < n; ++i) {
        if ((r = log_new(&map[i])))
            goto end;

        map[i].dup = false;
        map[i].srt = true;
        map[i].cmp.run = iasc;
    }

    for (int i = 0; i < sim->msh->hxd.len; ++i) {
        int *vtx = sim->msh->hxd.dat[i].vtx;

        for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 8; ++k)
                if (vtx[k] < vtx[j] && !log_add(&map[vtx[j]], vtx[k]))
                    z += 1;
    }

    if ((r = mtx_new(&fem->prv.ell, ((struct smtx_pps){n, z}))))
        goto end;

    for (int i = 0, e = 0; i < n; ++i) {
        fem->prv.ell.ia[i] = e;

        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e++)
            fem->prv.ell.ja[e] = j;
    }

    fem->prv.ell.ia[n] = z;

    if ((r = vec_new(&fem->prv.vec, n)))
        goto end;

    switch (sim->mod) {
        case SIM_HYP:
            if ((r = mtx_new(&fem->prv.hyp, fem->prv.ell.pps)))
                goto end;

            memcpy(fem->prv.hyp.ia, fem->prv.ell.ia, sizeof(int) * (n + 1));
            memcpy(fem->prv.hyp.ja, fem->prv.ell.ja, sizeof(int) * z);

            if ((r = mtx_new(&fem->prv.pbc, fem->prv.ell.pps)))
                goto end;

            memcpy(fem->prv.pbc.ia, fem->prv.ell.ia, sizeof(int) * (n + 1));
            memcpy(fem->prv.pbc.ja, fem->prv.ell.ja, sizeof(int) * z);

            break;
        case SIM_PBC:
            if ((r = mtx_new(&fem->prv.pbc, fem->prv.ell.pps)))
                goto end;

            memcpy(fem->prv.pbc.ia, fem->prv.ell.ia, sizeof(int) * (n + 1));
            memcpy(fem->prv.pbc.ja, fem->prv.ell.ja, sizeof(int) * z);

            break;
        case SIM_ELL:
            break;
    }

end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}
