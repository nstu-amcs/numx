#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <stdbool.h>

#include <numx/pde/slv.h>
#include <numx/vec/iss.h>

/** Finite Element Method simulation solver. */
struct fem
{
    struct slv slv;

    struct fem_ops
    {
        /** Solution mode. */
        enum
        {
            FEM_STD, // standard
            FEM_NON, // nonlinear
            FEM_HMC, // harmonic
        } mod;

        /** Basis functions. */
        enum
        {
            FEM_BSS_LIN, // triliniear
        } bss;

        /** Time discretization strategy. */
        enum
        {
            FEM_TDD_I2S, // implicit 2-layered
            FEM_TDD_I3S, // implicit 3-layered
            FEM_TDD_I4S, // implicit 4-layered
        } tdd;

        /** Options for nonlinear system solver. */
        struct
        {
            /**
             *  Nonlinear iteration callback (user defined).
             *
             *  Called for each nonlinear iteration.
             *  Doesn't called at all for standard and
             *  harmonic solution modes.
             *
             *  @param sim - simulation
             */
            struct
            {
                void *ctx;
                void (*run)(void *ctx, struct sim *sim);
            } itr;

            int    max;
            double err;
            double rlx;

            /** Runtime data made available by solver. */
            struct
            {
                int    itr;
                double err;
            } rt;
        } non;

        /** Options for linear system solver. */
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

    /** Private data for internal usage. */
    struct
    {
        struct smtx ell;
        struct smtx pbc;
        struct smtx hyp;
        struct vec  vec;
    } prv;
};

int fem_new(struct fem *fem);
int fem_cls(struct fem *fem);

int fem_exe(void *ctx, struct sim *sim);

#endif // NUMX_PDE_FEM_H
