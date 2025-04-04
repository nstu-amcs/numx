#ifndef NUMX_PDE_FEM_PRV_H
#define NUMX_PDE_FEM_PRV_H

#include <numx/pde/sim.h>
#include <numx/vec/mtx.h>

int fem_slv(struct sim *sim);

int fem_ell_slv(struct sim *sim);
int fem_ell_lin_std_slv(struct sim *sim, struct vec *wgt);
int fem_ell_lin_non_slv(struct sim *sim, struct vec *wgt);
int fem_ell_lin_hmc_slv(struct sim *sim, struct vec *wgt);

#endif // NUMX_PDE_FEM_PRV_H
