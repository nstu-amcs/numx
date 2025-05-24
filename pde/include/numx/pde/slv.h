#ifndef NUMX_PDE_SLV_H
#define NUMX_PDE_SLV_H

#include <numx/pde/msh.h>
#include <numx/vec/iss.h>
#include <numx/vec/vec.h>

#include <stdint.h>

#define NON_LAM(N) (N & (0b10000000))
#define NON_GAM(N) (N & (0b01000000))
#define NON_BET(N) (N & (0b00100000))
#define NON_SIG(N) (N & (0b00010000))
#define NON_CHI(N) (N & (0b00001000))
#define NON_SRC(N) (N & (0b00000100))
#define NON_TTA(N) (N & (0b00000010))
#define NON_EXT(N) (N & (0b00000001))

#define NON_ELL_MTX(N) ((N) & 0b11100000)
#define NON_ELL_VEC(N) ((N) & 0b00100111)
#define NON_ELL_ALL(N) (((N) & 0b11100000) && ((N) & 0b00100111))

#define NON_LAM_SET(N) (N |= (0b10000000))
#define NON_GAM_SET(N) (N |= (0b01000000))
#define NON_BET_SET(N) (N |= (0b00100000))
#define NON_SIG_SET(N) (N |= (0b00010000))
#define NON_CHI_SET(N) (N |= (0b00001000))
#define NON_SRC_SET(N) (N |= (0b00000100))
#define NON_TTA_SET(N) (N |= (0b00000010))
#define NON_EXT_SET(N) (N |= (0b00000001))

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
    struct slv_ops
    {
        /** Time discretization strategy. */
        enum tdd_mod
        {
            TDD_I2S = 2, // implicit 2-layered
            TDD_I3S = 3, // implicit 3-layered
            TDD_I4S = 4, // implicit 4-layered
        } tdd;

        /** Options for nonlinear system solver. */
        struct
        {
            enum non_mod
            {
                NON_FPI, // fixed-point iteration
                NON_NEW, // Newton's linearization
            } mod;

            uint8_t map;

            struct non_ops
            {
                int    max; // maximum number of iterations
                double err; // convergence tolerance
                bool   rlx; // enable relaxation

                /**
                 *  Nonlinear iteration callback (user defined).
                 *
                 *  Called for each nonlinear iteration.
                 */
                struct
                {
                    void *ctx;
                    void (*run)(void *ctx, struct non_ops *ops);
                } itr;

                enum
                {
                    DIF_NUM,
                    DIF_GIV,
                } dif;

                /** Runtime data made available by solver. */
                struct
                {
                    int    itr; // current iteration
                    double err; // current error
                    double rlx; // optimal relaxation factor
                } run;
            } ops;
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

        int    bs; // number of buffered layers (1 - 4)
        int    ti; // current time iteration
        double tv; // current time value
    } run;
} slv;

#endif // NUMX_PDE_SLV_H
