#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/pde/geo.h>

typedef double (*fun)(double, struct vtx *);

stdx_def_cut(fun_cut, fun);

/** Boundary condition. */
struct cnd
{
    enum cnd_mod
    {
        CND_DIR, // Dirichlet (1)
        CND_NEU, // Neumann (2)
        CND_ROB  // Robin (3)
    } mod;

    union
    {
        struct
        {
            fun tmp;
        } dir;

        struct
        {
            fun tta;
        } neu;

        struct
        {
            fun tmp;
            fun bet;
        } rob;
    } pps;
};

/** Read condition from given buffer using provided dataset. */
int cnd_get(struct cnd *cnd, const char *buf, struct fun_cut *dat);

#endif // NUMX_PDE_CND_H
