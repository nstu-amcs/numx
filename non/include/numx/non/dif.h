#ifndef NUMX_NON_DIF_H
#define NUMX_NON_DIF_H

#include <numx/vec/vec.h>

double pdif(void *ctx, double (*fun)(void *, struct vec *), int var, double hop, struct vec *vtx);

#endif // NUMX_NON_DIF_H
