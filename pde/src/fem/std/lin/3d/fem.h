#ifndef NUMX_PDE_PRV_FEM_STD_LIN_H
#define NUMX_PDE_PRV_FEM_STD_LIN_H

#include "../../ctx.h"

int fem_std_lin_slv(struct sim *sim, struct fem_std_ctx *ctx);

int fem_std_lin_asm(struct sim *sim, struct fem_std_ctx *ctx);
int fem_std_lin_new(struct sim *sim, struct fem_std_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_STD_LIN_H
