#ifndef NUMX_PDE_PRV_FEM_HMC_H
#define NUMX_PDE_PRV_FEM_HMC_H

#include <numx/pde/sim.h>

struct fem_hmc_ctx
{
    struct smtx mtx;
    struct vec  vec;
    struct vec  wgt;
};

#endif // NUMX_PDE_PRV_FEM_HMC_H
