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

int fem_lin_ell_slv(struct sim *sim);
int fem_lin_pbc_slv(struct sim *sim);
int fem_lin_hyp_slv(struct sim *sim);

int fem_lin_ell_asm(struct sim *sim);

int fem_lin_mtx_lam_asm(struct sim *sim, struct smtx *mtx);
int fem_lin_mtx_gam_asm(struct sim *sim, struct smtx *mtx);
int fem_lin_mtx_rob_asm(struct sim *sim, struct smtx *mtx);
int fem_lin_mtx_sig_asm(struct sim *sim, struct smtx *mtx);
int fem_lin_mtx_chi_asm(struct sim *sim, struct smtx *mtx);
int fem_lin_mtx_dir_asm(struct sim *sim, struct smtx *mtx);

int fem_lin_vec_src_asm(struct sim *sim, struct vec *vec);
int fem_lin_vec_neu_asm(struct sim *sim, struct vec *vec);
int fem_lin_vec_rob_asm(struct sim *sim, struct vec *vec);
int fem_lin_vec_dir_asm(struct sim *sim, struct vec *vec);

#endif // NUMX_PDE_FEM_PRV_H
