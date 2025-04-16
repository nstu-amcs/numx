#ifndef NUMX_PDE_SLV_H
#define NUMX_PDE_SLV_H

#include <numx/pde/msh.h>
#include <numx/vec/vec.h>

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
} slv;

#endif // NUMX_PDE_SLV_H
