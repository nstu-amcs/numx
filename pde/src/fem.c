#include <assert.h>
#include <numx/pde/sim.h>

int fem_slv(struct sim *sim)
{
    assert(sim);

    switch (sim->mod) {
        case SIM_ELL:
            return fem_ell_slv(sim);
            break;
        case SIM_PBC:
            return fem_pbc_slv(sim);
            break;
        case SIM_HYP:
            return fem_hyp_slv(sim);
            break;
    }

    return 0;
}
