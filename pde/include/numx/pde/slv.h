#ifndef NUMX_PDE_SLV_H
#define NUMX_PDE_SLV_H

#include <numx/msh/umsh.h>
#include <numx/vec/iss.h>
#include <numx/vec/vec.h>

#include <stdint.h>

struct sim;

/**
 * @brief Context passed to solution approximation function.
 */
struct apx_fun_ctx
{
    struct sim *sim; // Simulation.
    struct vec *wgt; // Solution (null for runtime).

    int var; // Variable to differentiate.
    int vtx; // Hinted vertex.
    int seg; // Hinted segment.
    int qud; // Hinted quadrangle.
    int hxd; // Hinted hexahedron.
};

/**
 * @brief Solution approximation function.
 *
 * @param ctx - context (`apx_fun_ctx`)
 * @param vtx - target point
 */
typedef double (*apx_fun)(void *ctx, struct vec *vtx);

/**
 * @brief Simulation solver.
 *
 * Each specific solver (such as FEM) will inherit general solver properties.
 */
typedef struct slv
{
    struct slv_ops
    {
        /** Time-domain discretization strategy. */
        enum tdd_mod
        {
            TDD_I2S = 2, // implicit 2-layered
            TDD_I3S = 3, // implicit 3-layered
            TDD_I4S = 4, // implicit 4-layered
        } tdd;

        /**
         * @brief Options for nonlinear system solver.
         */
        struct
        {
            bool enable; // Enable nonlinear mode.

            struct non_ops
            {
                int    max; // Maximum number of iterations.
                double err; // Convergence tolerance.

                bool rlx; // Enable relaxation.
                bool new; // Enable Newton's linearization tweaks.

                /**
                 *  @brief Initial field approximation.
                 *
                 *  If set, called once during startup.
                 *  If not set, initial field is set to zero.
                 */
                struct
                {
                    void *ctx;
                    void (*run)(void *ctx, struct vec *wgt);
                } ini_cbk;

                /**
                 *  @brief Nonlinear iteration callback (user defined).
                 *
                 *  Called for each nonlinear iteration.
                 */
                struct
                {
                    void *ctx;
                    void (*run)(void *ctx, struct non_ops *ops);
                } itr_cbk;

                enum
                {
                    DIF_NUM, // Differentiate numerically.
                    DIF_GIV, // Differentiate using given functions.
                } dif;       // Differentiation options.

                /** Runtime data made available by solver. */
                struct
                {
                    int    itr; // Current nonlinear iteration.
                    double err; // Current nonlinear error.
                    double rlx; // Optimal nonlinear relaxation factor.
                } run;
            } ops;
        } non;

        /**
         * @brief Options for linear system solver.
         */
        struct
        {
            enum iss_mod mod; // Linear system solver.

            union
            {
                struct iss_jac_ops jac;
                struct iss_rlx_ops rlx;
                struct iss_gmr_ops gmr;
                struct iss_bcg_ops bcg;
            } ops; // Solver-specific options.
        } iss;
    } ops;

    /**
     * @brief Execute the solver (implementation defined).
     */
    int (*exe)(struct sim *sim);

    /**
     *  @brief Solution approximation utilities (implementation defined).
     */
    struct
    {
        // Approximate the solution at the given point.
        apx_fun run;

        // Approximate the solution partial derivative at the given point.
        // Respected variable is passed via context `var` field.
        apx_fun dif;
    } apx;

    /**
     *  @brief Solution callback (user defined).
     *
     *  Called for each time layer. Called once
     *  for elliptic equations.
     */
    struct
    {
        void *ctx;
        void (*run)(void *ctx, struct sim *sim);
    } itr_cbk;

    /**
     *  @brief Runtime data made available by the solver.
     *
     *  Updated on each time layer and not
     *  available after solver exits.
     */
    struct
    {
        struct vec *wgt[4]; // buffered solution

        int    bs; // number of buffered layers (1 - 4)
        int    ti; // current time iteration
        double tv; // current time value
    } run;
} slv;

#endif // NUMX_PDE_SLV_H
