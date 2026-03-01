#ifndef NUMX_NON_DIF_H
#define NUMX_NON_DIF_H

#include <numx/non/fun.h>
#include <numx/vec/vec.h>

#define NUMX_NON_DIF_HOP 0.001

/**
 * @brief Partial differentiation options.
 */
typedef struct dif_ops
{
    struct vec *vtx; // fixed variables
    int         var; // respected variable

    double hop; // approximation step

    /**
     *  Tweak function.
     *
     *  Must increase respected variable by given hop.
     *  Can be used to differentite with respect to context.
     */
    void (*twk)(void *ctx, double hop, struct dif_ops *ops);
} dif_ops;

/**
 *  Three-point midpoint partial differentiation.
 */
double dif_tpm(void *ctx, mfun fun, struct dif_ops *ops);

#endif // NUMX_NON_DIF_H
