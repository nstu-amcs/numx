#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <numx/pde/sim.h>

struct fem_ops
{
    enum
    {
        FEM_STD,
        FEM_NON,
        FEM_HMC,
    } mod;

    enum
    {
        FEM_LIN,
    } bss;
};

struct fem_ell_ops
{
    struct fem_ops ops;
};

struct fem_pbc_ops
{
    struct fem_ops ops;
};

struct fem_hyp_ops
{
    struct fem_ops ops;
};

int fem_ell_slv(struct sim *sim, struct vec *q, struct fem_ell_ops ops);
int fem_pbc_slv(struct sim *sim, struct vec *q, struct fem_pbc_ops ops);
int fem_hyp_slv(struct sim *sim, struct vec *q, struct fem_hyp_ops ops);

int fem_slv(struct sim *sim, struct vec *q, struct fem_ops ops);

#endif // NUMX_PDE_FEM_H
