#ifndef NUMX_NON_NON_H
#define NUMX_NON_NON_H

#include <numx/com/cut.h>
#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>

enum non_new_mod
{
    EXC,
    CON
};

struct non_new_itr
{
    int k;

    double del;
    double err;

    struct vec *x;
};

struct non_new_ops
{
    enum non_new_mod mod;

    int    hem;
    double eps;
    double hop;

    struct non_new_itr *itr;

    struct vcap *cbk;
    struct jmtx *jac;
};

int non_new_slv(struct pcut *f, struct vec *x, struct non_new_ops o);

#endif // NUMX_NON_NON_H
