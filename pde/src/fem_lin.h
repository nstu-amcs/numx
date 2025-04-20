#ifndef NUMX_PDE_FEM_PRV_H
#define NUMX_PDE_FEM_PRV_H

#include <numx/pde/sim.h>
#include <numx/vec/mtx.h>

extern const double C;

extern const double G[2][2];
extern const double M[2][2];

extern const double GN[2][2];
extern const double MN[2][2][2];

extern const int MU[8];
extern const int NU[8];
extern const int TT[8];

int fem_lin_slv(struct sim *sim);
int fem_lin_asm(struct sim *sim);

#endif // NUMX_PDE_FEM_PRV_H
