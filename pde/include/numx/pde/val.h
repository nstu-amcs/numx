#ifndef NUMX_PDE_VAL_H
#define NUMX_PDE_VAL_H

#include <numx/com/cut.h>

struct sim;

typedef double (*fun)(struct sim *s, int vtx);

typedef struct val
{
    enum
    {
        VAL_NUM,
        VAL_FUN,
    } type;

    union
    {
        double num;
        fun    fun;
    } as;
} val;

cut_def(val_cut, val);

#endif // NUMX_PDE_VAL_H
