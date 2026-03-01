#ifndef NUMX_PDE_FEM_H
#define NUMX_PDE_FEM_H

#include <numx/pde/slv.h>

/**
 * @brief Finite Element Method simulation solver.
 */
struct fem
{
    /**
     * @brief Base solver properties.
     */
    struct slv slv;

    /**
     * @brief Implementation specific properties.
     */
    struct fem_ops
    {
        /**
         * @brief Basis functional space.
         */
        enum
        {
            FEM_BFS_LIN, // liniear
        } bfs;
    } ops;
};

/**
 * @brief Initialize the solver with defaults.
 */
int fem_new(struct fem *fem);

/**
 * @brief Execute the solver (implemetation entrypoint).
 */
int fem_exe(struct sim *sim);

#endif // NUMX_PDE_FEM_H
