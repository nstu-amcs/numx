#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <numx/pde/slv.h>

/** Finite Element Method simulation solver. */
struct fem
{
    struct slv slv;

    struct fem_ops
    {
        /** Basis functions. */
        enum
        {
            FEM_BSS_LIN, // triliniear
        } bss;
    } ops;
};

int fem_new(struct fem *fem);
int fem_exe(struct sim *sim);

#endif // NUMX_PDE_FEM_H
