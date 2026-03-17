/**
 * Approximation and interpolation subroutines.
 */

#ifndef NUMX_NON_APX_H
#define NUMX_NON_APX_H

#include <numx/com/cut.h>
#include <numx/vec/mtx.h>

/**
 *  @brief Interpolation function context.
 */
struct int_fun_ctx
{
    struct vec  *x; // Interpolation points.
    struct imtx *k; // Coefficient matrix.

    int prv; // Last accessed interval index.
};

/**
 *  @brief Cubic spline interpolation.
 *
 *  Coefficient matrix stores the result in the following form:
 *
 *  ```
 *   --                   --
 *   | a0 a1 a2 ... a(n-1) |
 *   | b0 b1 b2 ... b(n-1) |
 *   | c0 c1 c2 ... c(n-1) |
 *   | d0 d1 d2 ... d(n-1) |
 *   --                   --
 *  ```
 *
 *  Then, function value can be approximated as:
 *
 *  $$
 *  S_j(x) = a_j + b_j(x - x_j) + c_j(x - x_j)^2 + d_j(x - x_j)^3
 *  $$
 *
 *  @param v interpolation points
 *  @param k coefficient matrix
 */
int apx_int_cub(struct vec *x, struct imtx *k);

/**
 *  @brief Cubic spline interpolation function.
 *
 *  This is an `mfun` wrapper over spline interpolation
 *  function.
 *
 *  @param ctx interpolation context (`int_fun_ctx`)
 *  @param x point to approximate
 */
double int_cub_fun(void *ctx, struct vec *x);

/**
 *  @brief Cubic spline differentiation function.
 *
 *  This is an `mdif` wrapper over splite interpolation
 *  function differentiation.
 *
 *  @param ctx interpolation context (`int_fun_ctx`)
 *  @param fun cubic spline function (`int_cub_fun`)
 *  @param vtx point to differentiate at
 *  @param ops differentiation opsions (unused)
 */
double int_cub_dif(void *ctx, mfun fun, struct vec *vtx, struct dif_ops *ops);

#endif // NUMX_NON_APX_H
