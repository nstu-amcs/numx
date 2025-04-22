#ifndef NUMX_PDE_SLV_H
#define NUMX_PDE_SLV_H

#include <numx/pde/msh.h>
#include <numx/vec/vec.h>
#include <numx/vec/iss.h>

struct sim;

struct apx_fun_ctx
{
    struct sim *sim; // simulation
    struct vec *wgt; // solution (null for runtime)

    int vtx; // hinted vertex
    int qud; // hinted quadrangle
    int hxd; // hinted hexahedron
};

/** Simulation solver. */
typedef struct slv
{
    struct
    {
        /** Time discretization strategy. */
        enum
        {
            TDD_I2S = 2, // implicit 2-layered
            TDD_I3S = 3, // implicit 3-layered
            TDD_I4S = 4, // implicit 4-layered
        } tdd;

        /** Options for nonlinear system solver. */
        struct
        {
            enum
            {
                NON_FPI, // fixed-point iteration
                NON_NEW, // Newton's linearization
            } mod;

            struct
            {
                bool fd; // field dependence

                /**
                 *  Nonlinear iteration callback (user defined).
                 *
                 *  Called for each nonlinear iteration.
                 *
                 *  @param sim - simulation
                 */
                struct
                {
                    void *ctx;
                    void (*run)(void *ctx, struct sim *sim);
                } itr;

                int    max; // maximum number of iterations
                double err; // convergence tolerance
                double rlx; // relaxation factor
            } ops;

            /** Runtime data made available by solver. */
            struct
            {
                int    itr; // current iteration
                double err; // current error
            } run;
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

    /**
     *  Execute solver (implementation defined).
     *
     *  @param sim - simulation
     */
    int (*exe)(struct sim *sim);

    /**
     *  Approximate solution at given point (implementation defined).
     *
     *  @param ctx - context
     *  @param vtx - target point
     */
    double (*apx)(struct apx_fun_ctx *ctx, struct vec *vtx);

    /**
     *  Solution callback (user defined).
     *
     *  Called for each time layer. Called once
     *  for elliptic equations.
     *
     *  @param sim - simulation
     */
    struct
    {
        void *ctx;
        void (*run)(void *ctx, struct sim *sim);
    } itr;

    /**
     *  Runtime data made available by solver.
     *
     *  Updated on each time layer and not
     *  available after solver exits.
     */
    struct
    {
        struct vec *wgt[4]; // buffered solution

        int tb; // number of buffered layers (1 - 4)
        int ti; // current time iteration
        double tv; // current time value
    } run;
} slv;

#endif // NUMX_PDE_SLV_H
