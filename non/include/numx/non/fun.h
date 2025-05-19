#ifndef NUMX_NON_FUN_H
#define NUMX_NON_FUN_H

#include <numx/vec/vec.h>

typedef double (*ufun)(void *, double);
typedef double (*mfun)(void *, struct vec *);

#endif // NUMX_NON_FUN_H
