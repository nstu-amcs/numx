#include <assert.h>
#include <numx/pde/fem.h>

int fem_slv(struct sim *sim, struct vec *q)
{
    assert(sim);
    assert(q);

    switch (sim->mod) {
        case SIM_ELL:
            return fem_ell_slv(sim, q);
            break;
        case SIM_PBC:
            return fem_pbc_slv(sim, q);
            break;
        case SIM_HYP:
            return fem_hyp_slv(sim, q);
            break;
    }

    return 0;
}
