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

            char dir[128]; // Export directory.
            char pfx[64];  // Export prefix.
            char sol[64];  // Primary solution name.

            /**
             *  @brief Export commons (defined by `mod`).
             *
             *  @param sim - simulation
             */
            int (*ini)(struct sim *sim);

            /**
             *  @brief Export vertex-wise solution (defined by `mod`).
             *
             *  @param sim simulation
             *  @param name solution name
             *  @param i solution index (0 for primary)
             *  @param s solution data
             */
            int (*put_v)(struct sim *sim, const char *name, int i, struct vec *s);

            /**
             *  @brief Export cell-wise solution (defined by `mod`).
             *
             *  @param sim simulation
             *  @param name solution name
             *  @param i solution index (0 for primary)
             *  @param s solution data
             */
            int (*put_c)(struct sim *sim, const char *name, int i, struct vec *s);
        } exp;

        /**
         * @brief Time discretization options.
         */
        struct
        {
            int num; // Number of time intervals.
            int ini; // Number of precomputed layers.

            double beg; // Initial time.
            double hop; // Time interval length.
        } tdd;

        /**
         * @brief Harmonic options.
         */
        struct
        {
            double frq; // Frequency.
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
 * @brief Export vertex-wise solution in CGNS format.
 */
int sim_exp_cgns_put_v(struct sim *sim, const char *name, int i, struct vec *s);

/**
 * @brief Export cell-wise solution in CGNS format.
 */
int sim_exp_cgns_put_c(struct sim *sim, const char *name, int i, struct vec *s);

/**
 * @brief Start the simulation.
 */
int sim_run(struct sim *sim);

#endif // NUMX_PDE_SIM_H
