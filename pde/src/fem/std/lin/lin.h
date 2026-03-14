#ifndef NUMX_PDE_PRV_FEM_STD_LIN_H
#define NUMX_PDE_PRV_FEM_STD_LIN_H

#include "../std.h"

double fem_std_lin_c2d_apx(void *ctx, struct vec *vtx);
double fem_std_lin_c2d_dif(void *ctx, struct vec *vtx);

double fem_std_lin_c3d_apx(void *ctx, struct vec *vtx);
double fem_std_lin_c3d_dif(void *ctx, struct vec *vtx);

int fem_std_lin_c2d_asm(struct sim *sim, struct fem_std_ctx *ctx);
int fem_std_lin_c3d_asm(struct sim *sim, struct fem_std_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_STD_LIN_H
