#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/com/cut.h>
#include <numx/pde/msh.h>
#include <numx/vec/vec.h>

struct sim;

struct val_fun_ctx
{
    struct sim *sim;
    int         vtx;
};

typedef struct val
{
    enum
    {
        VAL_NUM, // constant
        VAL_FUN, // function
    } type;

    union
    {
        double num;
        double (*fun)(void *ctx, struct vec *vtx); // function of space, time and field
    } as;

    struct
    {
        double (*dif)(void *ctx, struct vec *vtx); // partial derivative with respect to field
    } ops;
} val;

/** Boundary condition. */
typedef struct cnd_bnd
{
    /** Boundary condition type. */
    enum cnd_bnd_type
    {
        CND_BND_DIR, // Dirichlet
        CND_BND_NEU, // Neumann
        CND_BND_ROB, // Robin
    } type;

    union
    {
        struct
        {
            val tgt; // field
        } dir;

        struct
        {
            val tta; // field flux
        } neu;

        struct
        {
            val bet; // robin coefficient
            val ext; // external field
        } rob;
    } pps;
} cnd_bnd;

/** Initial condition. */
typedef struct cnd_ini
{
    val tgt; // field
} cnd_ini;

cut_def(val_cut, val);
cut_def(cnd_bnd_cut, cnd_bnd);
cut_def(cnd_ini_cut, cnd_ini);

#endif // NUMX_PDE_CND_H
