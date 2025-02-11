#ifndef NUMX_NON_DIF_H
#define NUMX_NON_DIF_H

#include <numx/vec/vec.h>

int pdif(double (*f)(struct vec*), int i, double hop, struct vec* p, double* r);

#endif  // NUMX_NON_DIF_H
