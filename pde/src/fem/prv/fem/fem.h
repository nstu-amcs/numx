#ifndef NUMX_PDE_PRV_FEM_H
#define NUMX_PDE_PRV_FEM_H

#include <numx/pde/sim.h>
#include <numx/vec/mtx.h>

#include <stdint.h>

struct fem_ctx
{
    struct smtx mtx;
    struct vec  vec;

    struct smtx sig;
    struct smtx chi;

    struct vec w0;
    struct vec w1;
    struct vec w2;
    struct vec w3;
};

#endif // NUMX_PDE_PRV_FEM_H
