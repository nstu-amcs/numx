#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <numx/pde/sim.h>

int fem_ell_slv(struct sim *sim, struct vec *q);
int fem_pbc_slv(struct sim *sim, struct vec *q);
int fem_hyp_slv(struct sim *sim, struct vec *q);

int fem_slv(struct sim *sim, struct vec *q);

#endif // NUMX_PDE_FEM_H
