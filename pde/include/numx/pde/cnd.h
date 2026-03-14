#ifndef NUMX_PDE_CND_H
#define NUMX_PDE_CND_H

#include <numx/com/cut.h>
#include <numx/non/fun.h>
#include <numx/vec/vec.h>

typedef struct sim sim;

/**
 * @brief Context passed to parameters given as functions.
 */
struct sim_fun_ctx
{
    struct sim *sim;
    void       *ctx; // User-defined context.

    int vtx; // Hinted vertex.
    int seg; // Hinted segment.
    int qud; // Hinted quadrangle.
    int hxd; // Hinted hexahedron.
};

/**
 * @brief Generic property value.
 */
typedef struct val
{
    enum
    {
        VAL_NUM, // Constant.
        VAL_FUN, // Function.
        VAL_HMC, // Harmonic.
    } type;

    union
    {
        double num; // Constant value.

        struct
        {
            void *ctx; // User-defined context.

            // Function of space, time and field.
            // The `sim_fun_ctx` will be passed as context with
            // user-defined context inside.
            mfun run;
        } fun;

        struct
        {
            mfun sin;
            mfun cos;
        } hmc;
    } as;

    struct
    {
        bool dep; // field dependence
        mfun dif; // partial derivative with respect to field
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
