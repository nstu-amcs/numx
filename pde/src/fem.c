#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/fem.h>

#include "fem.h"

int fem_new(struct fem *fem)
{
    assert(fem);

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

    fem->ops.tdd = FEM_TDD_I4S;

    return 0;
}

int fem_cls(struct fem *fem)
{
    assert(fem);
    return 0;
}

static int fem_ini(struct sim *sim);

int fem_slv(struct sim *sim)
{
    assert(sim);

    if (fem_ini(sim))
        return -1;

    switch (sim->mod) {
        case SIM_ELL:
            return fem_ell_slv(sim);
        case SIM_PBC:
            errno = ENOTSUP;
            return -1;
        case SIM_HYP:
            errno = ENOTSUP;
            return -1;
    }

    return 0;
}

static int fem_ini(struct sim *sim)
{
    int n = sim->msh->vtx.len;
    int z = 0;
    int r = 0;

    struct ilog *map = malloc(sizeof(struct ilog) * n);

    if (!map) {
        r = -1;
        goto ini_end;
    }

    for (int i = 0; i < n; ++i) {
        if ((r = log_new(&map[i])))
            goto ini_end;

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

    if ((r = mtx_new(&sim->fem->mtx, ((struct smtx_pps){n, z}))))
        goto ini_end;

    for (int i = 0, e = 0; i < n; ++i) {
        sim->fem->mtx.ia[i] = e;

        log_rst(&map[i]);

        for (int j = 0; !log_adv(&map[i], &j); e++)
            sim->fem->mtx.ja[e] = j;
    }

    sim->fem->mtx.ia[n] = z;

    if ((r = vec_new(&sim->fem->vec, n)))
        goto ini_end;

ini_end:
    for (int i = 0; i < n; ++i)
        log_cls(&map[i]);

    free(map);

    return r;
}
