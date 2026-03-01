#ifndef NUMX_NON_FUN_H
#define NUMX_NON_FUN_H

#include <numx/vec/vec.h>

/**
 * @brief Univariate function.
 */
typedef double (*ufun)(void *, double);

/**
 * @brief Multivariate function.
 */
typedef double (*mfun)(void *, struct vec *);

#endif // NUMX_NON_FUN_H
