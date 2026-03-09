/**
 * Approximation and interpolation subroutines.
 */

#ifndef NUMX_NON_APX_H
#define NUMX_NON_APX_H

#include <numx/com/cut.h>
#include <numx/vec/mtx.h>

/**
 * @brief Cubic Hermite spline interpolation.
 */
int apx_cub(struct vec *x, struct imtx *k);

#endif // NUMX_NON_APX_H
