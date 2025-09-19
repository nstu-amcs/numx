#ifndef NUMX_INV_INV_H
#define NUMX_INV_INV_H

#include <numx/vec/vec.h>

struct rec
{
    double wgt;
    double (*eps)(struct vec *);
};

#endif // NUMX_INV_INV_H
