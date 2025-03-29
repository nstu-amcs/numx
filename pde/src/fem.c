#include <assert.h>
#include <numx/pde/fem.h>

int fem_slv(struct sim *sim, struct vec *q, struct fem_ops ops)
{
    assert(sim);
    assert(q);

    switch (sim->mod) {
        case SIM_ELL:
            return fem_ell_slv(sim, q, *((struct fem_ell_ops *)&ops));
            break;
        case SIM_PBC:
            return fem_pbc_slv(sim, q, *((struct fem_pbc_ops *)&ops));
            break;
        case SIM_HYP:
            return fem_hyp_slv(sim, q, *((struct fem_hyp_ops *)&ops));
            break;
    }

    return 0;
}
