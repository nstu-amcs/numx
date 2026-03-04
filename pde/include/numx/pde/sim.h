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
    val sig; // capacity coefficient
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

/**
 * @brief Simulation.
 */
typedef struct sim
{
    /**
     * @brief Equation type.
     */
    enum
    {
        SIM_ELL, // elliptic
        SIM_PBC, // parabolic
        SIM_HYP, // hyperbolic
    } eqn;

    /**
     * @brief Solution mode.
     */
    enum
    {
        SIM_STD, // standard
        SIM_HMC, // harmonic
    } mod;

    /**
     * @brief Simulation options.
     */
    struct sim_ops
    {
        /**
         * @brief User-defined functions.
         */
        struct
        {
            char dir[128]; // usr directory
            char pfx[64];  // usr prefix

            void *hdl; // dynamic-library handler
        } usr;

        /**
         * @brief Export options.
         */
        struct
        {
            enum
            {
                SIM_EXP_CGNS, // CGNS
            } mod;

            char dir[128]; // export directory
            char pfx[64];  // export prefix

            /**
             *  @brief Export commons (simulation defined).
             *
             *  @param sim - simulation
             */
            int (*ini)(struct sim *sim);

            /**
             *  @brief Export runtime solution (simulation defined).
             *
             *  @param sim - simulation
             */
            int (*put)(struct sim *sim);
        } exp;

        /**
         * @brief Time discretization options.
         */
        struct
        {
            int num; // number of time intervals
            int ini; // number of precomputed layers

            double beg; // initial time
            double hop; // time interval length
        } tdd;

        /**
         * @brief Harmonic options.
         */
        struct
        {
            double frq; // frequency
        } hmc;
    } ops;

    struct umsh *msh; // active mesh
    struct slv  *slv; // active solver

    struct mat_cut mat; // materials
    struct val_cut src; // sources
    struct obj_cut obj; // objects
    struct bnd_cut bnd; // boundaries

    struct cnd_ini_cut cnd_ini; // initial conditions
    struct cnd_bnd_cut cnd_bnd; // boundary conditions
} sim;

/**
 * @brief Initialize simulation with defaults.
 *
 * Must be called before any other method.
 */
int sim_new(struct sim *sim);

/**
 * @brief Finish simulation and free memory.
 */
int sim_cls(struct sim *sim);

/**
 * @brief Import simulation from the CGNS file.
 */
int sim_imp_cgns(struct sim *sim, const char *gns);

/**
 * @brief Import simulation from the Elmer file.
 */
int sim_imp_elm(struct sim *sim, const char *sif);

/**
 * @brief Import simulation from the Telma file.
 */
int sim_imp_tel(struct sim *sim, const char *tel);

/**
 * @brief Export commons in CGNS format.
 */
int sim_exp_cgns_ini(struct sim *sim);

/**
 * @brief Export solution in CGNS format.
 */
int sim_exp_cgns_put(struct sim *sim);

/**
 * @brief Start the simulation.
 */
int sim_run(struct sim *sim);

#endif // NUMX_PDE_SIM_H
