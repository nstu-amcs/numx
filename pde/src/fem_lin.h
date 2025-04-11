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

double fem_lin_apx(void *ctx, struct sim *s, struct vec *wgt, int hxd, struct vtx *vtx);

int fem_lin_ell_mtx_asm(struct sim *sim, struct vec *wgt);
int fem_lin_ell_vec_asm(struct sim *sim, struct vec *wgt);
int fem_lin_pbc_mtx_asm(struct sim *sim, struct vec *wgt);
int fem_lin_hyp_mtx_asm(struct sim *sim, struct vec *wgt);

int fem_lin_ell_slv(struct sim *sim);
int fem_lin_pbc_slv(struct sim *sim);
int fem_lin_hyp_slv(struct sim *sim);

int asm_mov_mtx(struct smtx *m, int i, int j, double v);

#endif // NUMX_PDE_FEM_PRV_H
