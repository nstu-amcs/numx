#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/com/cut.h>
#include <numx/pde/msh.h>
#include <numx/pde/val.h>

typedef struct cnd_bnd
{
    enum cnd_bnd_type
    {
        CND_BND_DIR,
        CND_BND_NEU,
        CND_BND_ROB,
    } type;

    union
    {
        struct
        {
            val tgt;
        } dir;

        struct
        {
            val tta;
        } neu;

        struct
        {
            val tgt;
            val bet;
            val src;
        } rob;
    } pps;
} cnd_bnd;

typedef struct cnd_ini
{
    val tgt;
} cnd_ini;

cut_def(cnd_bnd_cut, cnd_bnd);
cut_def(cnd_ini_cut, cnd_ini);

#endif // NUMX_PDE_CND_H
