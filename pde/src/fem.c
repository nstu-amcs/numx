#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>

#include "fem_lin.h"

int fem_new(struct fem *fem)
{
    assert(fem);

    fem->slv.exe = fem_exe;

    fem->ops.mod = FEM_STD;
    fem->ops.bss = FEM_BSS_LIN;

    fem->ops.iss.mod = ISS_BCG;
    fem->ops.iss.ops.bcg.con.sm = 0;
    fem->ops.iss.ops.bcg.ops.err = 1e-10;
    fem->ops.iss.ops.bcg.ops.itr.run = 0;
    fem->ops.iss.ops.bcg.ops.max = 500;

    fem->ops.non.mod = NON_FPI;
    fem->ops.non.ops.itr.ctx = 0;
    fem->ops.non.ops.itr.run = 0;
    fem->ops.non.ops.max = 50;
    fem->ops.non.ops.err = 1e-10;
    fem->ops.non.ops.rlx = 1;

    fem->ops.tdd = FEM_TDD_I2S;

    return 0;
}

static int fem_ini(struct sim *sim);

int fem_exe(struct sim *sim)
{
    if (fem_ini(sim))
        return -1;

    if (sim->ops.exp.ini(sim))
        return -1;

    switch (((struct fem *)sim->slv)->ops.bss) {
        case FEM_BSS_LIN:
            return fem_lin_slv(sim);
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

    if ((r = mtx_new(&fem->prv.mtx, ((struct smtx_pps){n, z}))))
        goto end;

    for (int i = 0, e = 0; i < n; ++i) {
        fem->prv.mtx.ia[i] = e;

        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e++)
            fem->prv.mtx.ja[e] = j;
    }

    fem->prv.mtx.ia[n] = z;

    if ((r = vec_new(&fem->prv.vec, n)))
        goto end;

    switch (sim->mod) {
        case SIM_HYP:
            if ((r = mtx_new(&fem->prv.chi, fem->prv.mtx.pps)))
                goto end;

            memcpy(fem->prv.chi.ia, fem->prv.mtx.ia, sizeof(int) * (n + 1));
            memcpy(fem->prv.chi.ja, fem->prv.mtx.ja, sizeof(int) * z);

        [[fallthrough]];
        case SIM_PBC:
            if ((r = mtx_new(&fem->prv.sig, fem->prv.mtx.pps)))
                goto end;

            memcpy(fem->prv.sig.ia, fem->prv.mtx.ia, sizeof(int) * (n + 1));
            memcpy(fem->prv.sig.ja, fem->prv.mtx.ja, sizeof(int) * z);

        [[fallthrough]];
        case SIM_ELL:
            break;
    }

end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}
