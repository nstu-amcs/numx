#ifndef NUMX_PDE_SIM_H
#define NUMX_PDE_SIM_H

#include <numx/pde/cnd.h>
#include <numx/pde/fem.h>
#include <numx/vec/iss.h>

/**
 * Material.
 *
 * Specifies the physical properties of an object.
 */
typedef struct mat
{
    val lam; // diffusion coefficient
    val gam; // reaction coefficient
    val sig;
    val chi;
} mat;

/**
 *  Object.
 *
 *  Physical group for volume elements.
 */
typedef struct obj
{
    int mat; // material
    int ini; // initial condition
    int src; // source field
} obj;

/**
 *  Boundary.
 *
 *  Physical group for face elements.
 */
typedef struct bnd
{
    int cnd; // boundary condition
} bnd;

cut_def(mat_cut, mat);
cut_def(obj_cut, obj);
cut_def(bnd_cut, bnd);

/** Simulation. */
typedef struct sim
{
    /** Simulation mode (equation type). */
    enum
    {
        SIM_ELL, // elliptic
        SIM_PBC, // parabolic
        SIM_HYP, // hyperbolic
    } mod;

    struct sim_ops
    {
        void *usr; // user-defined functions

        /** Export options. */
        struct
        {
            enum
            {
                SIM_EXP_GNS, // cgns
            } mod;

            char dir[128]; // export directory
            char pfx[64];  // export prefix

            /**
             *  Export commons (simulation defined).
             *
             *  @param sim - simulation
             */
            struct
            {
                void *ctx;
                int (*run)(void *ctx, struct sim *sim);
            } ini;

            /**
             *  Export solution (simulation defined).
             *
             *  @param sim - simulation
             */
            struct
            {
                void *ctx;
                int (*run)(void *ctx, struct sim *sim);
            } put;
        } exp;

        /** Time discretization options. */
        struct
        {
            int num; // number of time intervals
            int hop; // time interval length
        } tdd;
    } ops;

    struct msh *msh; // active mesh
    struct slv *slv; // active solver

    struct mat_cut mat; // materials
    struct val_cut src; // sources
    struct obj_cut obj; // objects
    struct bnd_cut bnd; // boundaries

    struct cnd_ini_cut cnd_ini; // initial conditions
    struct cnd_bnd_cut cnd_bnd; // boundary conditions

    /**
     *  Runtime data made available by solver.
     *
     *  Updated on each time layer and not
     *  available after solver exits.
     */
    struct
    {
        int ti; // current time iteration
        int tv; // current time value

        struct vec *wgt; // current solution
    } rt;
} sim;

int sim_new(struct sim *sim);
int sim_cls(struct sim *sim);

/** Import simulation from CGNS file. */
int sim_imp_gns(struct sim *sim, const char *gns);

/** Import simulation from Elmer file. */
int sim_imp_elm(struct sim *sim, const char *sif);

/** Export commons in CGNS format. */
int sim_exp_ini_gns(void *ctx, struct sim *sim);

/** Export solution in CGNS format. */
int sim_exp_put_gns(void *ctx, struct sim *sim);

/** Start simulation. */
int sim_run(struct sim *sim);

#endif // NUMX_PDE_SIM_H
