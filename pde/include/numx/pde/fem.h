#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <stdbool.h>

#include <numx/vec/iss.h>

struct sim;

struct fem
{
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
            FEM_BSS_LIN,
        } bss;

        enum
        {
            FEM_TDD_I2S,
            FEM_TDD_I3S,
            FEM_TDD_I4S,
        } tdd;

        struct
        {
            int    max;
            double err;
            double rlx;
        } non;

        struct
        {
            enum iss_mod mod;
            union
            {
                struct iss_jac_ops jac;
                struct iss_rlx_ops rlx;
                struct iss_bcg_ops bcg;
            } ops;
        } iss;
    } ops;

    struct smtx mtx;
    struct vec  vec;
};

int fem_new(struct fem *fem);
int fem_cls(struct fem *fem);

#endif // NUMX_PDE_FEM_H
