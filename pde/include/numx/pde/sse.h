#ifndef NUMX_PDE_SSE_H
#define NUMX_PDE_SSE_H

#include <numx/pde/cnd.h>
#include <numx/pde/geo.h>
#include <numx/vec/iss.h>
#include <numx/vec/vec.h>
#include <stdx/cut.h>

/** Face's context. */
struct fctx
{
    struct cnd cnd;
};

int fctx_get(void *ctx, int n, ...);

/** Finite element's context. */
struct ectx
{
    fun lam;
    fun gam;
    fun ext;
};

int ectx_get(void *ctx, int n, ...);

/** Options for boundary problem solver. */
struct sse_ops
{
    int    max; // Maximum number of iterations (set 1 for linear problems)
    double eps;

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
};

/** Options for finite difference boundary problem solver. */
struct sse_fdm_ops
{
    struct sse_ops ops;
};

/** Options for finite element boundary problem solver. */
struct sse_fem_ops
{
    struct sse_ops ops;

    enum fem_mod
    {
        FEM_LIN // Trilinear basis
    } mod;
};

/// @brief Finite Difference boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fdm_slv(struct obj *obj, struct vec *x, struct sse_fdm_ops ops);

/// @brief Finite Element boundary problem solver.
/// @param o domain object
/// @param x solution vector
int pde_sse_fem_slv(struct obj *obj, struct vec *x, struct sse_fem_ops ops);

#endif // NUMX_PDE_SSE_H
