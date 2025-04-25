#ifndef NUMX_NON_OPM_H
#define NUMX_NON_OPM_H

#include <numx/vec/vec.h>

#define NUMX_NON_OPM_EPS 0.001

typedef struct opm_ops
{
    struct vec *vtx; // fixed variables
    int         var; // respected variable

    double eps; // precision
    double beg; // range start
    double end; // range end

    /**
     *  Tweak function.
     *
     *  Must replace respected variable with given value.
     *  Can be used to optimize with respect to context.
     * */
    void (*twk)(void *ctx, double val, struct opm_ops *ops);
} opm_ops;

/**
 *  Bisection univariable local optimization.
 */
double opm_loc_bis(void *ctx, double (*fun)(void *, struct vec *), struct opm_ops *ops);

#endif // NUMX_NON_OPM_H
