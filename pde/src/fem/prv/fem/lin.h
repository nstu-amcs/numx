#ifndef NUMX_PDE_PRV_FEM_LIN_H
#define NUMX_PDE_PRV_FEM_LIN_H

#include <numx/pde/sim.h>
#include <numx/vec/mtx.h>

#include <prv/fem/fem.h>

int fem_lin_slv(struct sim *sim, struct fem_ctx *ctx);
int fem_lin_asm(struct sim *sim, struct fem_ctx *ctx);
int fem_lin_new(struct sim *sim, struct fem_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_LIN_H
