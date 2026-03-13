/**
 * Approximation and interpolation subroutines.
 */

#ifndef NUMX_NON_APX_H
#define NUMX_NON_APX_H

#include <numx/com/cut.h>
#include <numx/vec/mtx.h>

/**
 * @brief Cubic Hermite spline interpolation.
 * 
 * Coefficient matrix stores the result in the following form:
 * 
 * ```
 *  --                   --
 *  | a0 a1 a2 ... a(n-1) |
 *  | b0 b1 b2 ... b(n-1) |
 *  | c0 c1 c2 ... c(n-1) |
 *  | d0 d1 d2 ... d(n-1) |
 *  --                   --
 * ```
 * 
 * @param v interpolation points
 * @param k coefficient matrix
 */
int apx_int_cub(struct vec *x, struct imtx *k);

#endif // NUMX_NON_APX_H
