#ifndef NUMX_PDE_SLV_H
#define NUMX_PDE_SLV_H

#include <numx/pde/msh.h>
#include <numx/vec/vec.h>

struct sim;

/** Simulation solver. */
typedef struct slv
{
    /**
     *  Execute solver (implementation defined).
     *
     *  @param sim - simulation
     */
    struct
    {
        void *ctx;
        int (*run)(void *ctx, struct sim *sim);
    } exe;

    /**
     *  Approximate solution at given point (implementation defined).
     *
     *  @param sim - simulation
     *  @param wgt - solution
     *  @param hxd - target element
     *  @param vtx - target point
     */
    struct
    {
        void *ctx;
        double (*run)(void *ctx, struct sim *s, struct vec *wgt, int hxd, struct vec *vtx);
    } apx;

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
