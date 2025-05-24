#ifndef NUMX_PDE_PRV_FEM_HMC_LIN_H
#define NUMX_PDE_PRV_FEM_HMC_LIN_H

#include "../ctx.h"

int fem_hmc_lin_slv(struct sim *sim, struct fem_hmc_ctx *ctx);
int fem_hmc_lin_asm(struct sim *sim, struct fem_hmc_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_HMC_LIN_H
