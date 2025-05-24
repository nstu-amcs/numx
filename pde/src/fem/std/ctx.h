#ifndef NUMX_PDE_PRV_FEM_STD_CTX_H
#define NUMX_PDE_PRV_FEM_STD_CTX_H

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

#endif // NUMX_PDE_PRV_FEM_STD_CTX_H
