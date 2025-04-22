#ifndef NUMX_NON_DIF_H
#define NUMX_NON_DIF_H

#include <numx/vec/vec.h>

struct dif_ops
{
    double hop; // approximation step

    /** Tweak function. */
    void (*twk)(void*, struct vec *, int, double);
};

double pdif(void *ctx, double (*fun)(void *, struct vec *), int var, struct vec *vtx, struct dif_ops *ops);

#endif // NUMX_NON_DIF_H
