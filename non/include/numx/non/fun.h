#ifndef NUMX_NON_FUN_H
#define NUMX_NON_FUN_H

#include <numx/vec/vec.h>

#define NUMX_NON_DIF_HOP 0.001

/**
 * @brief Univariate function.
 */
typedef double (*ufun)(void *, double);

/**
 * @brief Multivariate function.
 */
typedef double (*mfun)(void *, struct vec *);

/**
 * @brief Partial differentiation options.
 */
typedef struct dif_ops
{
    int    var; // Respected variable.
    double hop; // Approximation step.

    /**
     *  Tweak function.
     *
     *  Must increase respected variable by given hop.
     *  Can be used to differentite with respect to context.
     */
    void (*twk)(void *ctx, double hop, struct dif_ops *ops);
} dif_ops;

/**
 *  @brief Partial differentiation function.
 *
 *  @param ctx underlying function context
 *  @param fun function to differentiate
 *  @param vtx point to differentiate at
 *  @param ops differentiation options
 */
typedef double (*mdif)(void *ctx, mfun fun, struct vec *vtx, struct dif_ops *ops);

/**
 *  Three-point midpoint partial numeric differentiation.
 */
double dif_tpm(void *ctx, mfun fun, struct vec *vtx, struct dif_ops *ops);

#endif // NUMX_NON_FUN_H
