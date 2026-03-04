#ifndef NUMX_PDE_PRV_FEM_HMC_H
#define NUMX_PDE_PRV_FEM_HMC_H

#include <numx/pde/sim.h>

struct fem_hmc_ctx
{
    struct smtx mtx;
    struct vec  vec;
    struct vec  wgt;
};

int fem_hmc_lin_slv(struct sim *sim, struct fem_hmc_ctx *ctx);
int fem_hmc_lin_asm(struct sim *sim, struct fem_hmc_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_HMC_H
