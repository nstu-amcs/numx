#ifndef NUMX_PDE_PRV_FEM_STD_H
#define NUMX_PDE_PRV_FEM_STD_H

#include <numx/pde/sim.h>
#include <numx/vec/mtx.h>

struct fem_std_ctx
{
    struct smtx mtx;
    struct smtx sig;
    struct smtx chi;
    struct vec  vec;

    struct vec w0;
    struct vec w1;
    struct vec w2;
    struct vec w3;

    struct vec tmp;
};

typedef int (*fem_asm_fun)(struct sim *, struct fem_std_ctx *);

double fem_std_lin_apx(struct apx_fun_ctx *ctx, vtx_ptr vtx);
int    fem_std_lin_asm(struct sim *sim, struct fem_std_ctx *ctx);

#endif // NUMX_PDE_PRV_FEM_STD_H
