#include <assert.h>

#include <numx/com/cmp.h>
#include <numx/com/log.h>
#include <numx/pde/sim.h>

#include "fem.h"

int fem_ell_slv(struct sim *sim)
{
    assert(sim);

    int        r = 0;
    struct vec wgt;

    if ((r = vec_new(&wgt, sim->fem->vec.n)))
        goto end;

    switch (sim->fem->ops.bss) {
        case FEM_BSS_LIN:
            switch (sim->fem->ops.mod) {
                case FEM_STD:
                    if ((r = fem_ell_lin_std_slv(sim, &wgt)))
                        goto end;

                    break;
                case FEM_NON:
                    if ((r = fem_ell_lin_non_slv(sim, &wgt)))
                        goto end;

                    break;
                case FEM_HMC:
                    errno = ENOTSUP;
                    return -1;
            }

            break;
        default:
            errno = ENOTSUP;
            return -1;
    }

    sim->ops.exp.put.run(sim->ops.exp.put.ctx, 2, 0, &wgt);

end:
    vec_cls(&wgt);

    return r;
}
