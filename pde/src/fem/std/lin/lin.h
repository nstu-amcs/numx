#ifndef NUMX_PDE_PRV_FEM_STD_LIN_H
#define NUMX_PDE_PRV_FEM_STD_LIN_H

#include "../std.h"

double fem_std_lin_c2d_apx(struct apx_fun_ctx *ctx, vtx_ptr vtx);
double fem_std_lin_c3d_apx(struct apx_fun_ctx *ctx, vtx_ptr vtx);

int    fem_std_lin_c2d_asm(struct sim *sim, struct fem_std_ctx *ctx);
int    fem_std_lin_c3d_asm(struct sim *sim, struct fem_std_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_STD_LIN_H