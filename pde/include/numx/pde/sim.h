#ifndef NUMX_PDE_SIM_H
#define NUMX_PDE_SIM_H

#include <numx/pde/cnd.h>
#include <numx/vec/iss.h>

typedef struct mat
{
    val lam;
    val gam;
    val sig;
    val chi;
} mat;

typedef struct obj
{
    int mat;
    int ini;
    int ext;
} obj;

typedef struct bnd
{
    int cnd;
} bnd;

stdx_def_cut(mat_cut, mat);
stdx_def_cut(obj_cut, obj);
stdx_def_cut(bnd_cut, bnd);

struct sim_ops
{
    struct
    {
        int    max;
        double eps;
        double omg;
    } non;

    struct
    {
        enum iss_mod mod;
        union
        {
            struct iss_jac_ops jac;
            struct iss_rlx_ops rlx;
            struct iss_bcg_ops bcg;
        } ops;
    } iss;
};

struct sim_sse_ops
{
    struct sim_ops ops;
};

struct sim_tde_ops
{
    struct sim_ops ops;

    int num;
    int hop;
};

struct sim
{
    struct
    {
        struct
        {
            char dir[64];
            char pfx[64];
        } msh;

        struct
        {
            char dir[64];
            char pfx[64];
        } exp;

        void *usr;
    } pps;

    enum sim_mod
    {
        SIM_SSE,
        SIM_TDE,
    } mod;

    union
    {
        struct sim_sse_ops sse;
        struct sim_tde_ops tde;
    } ops;

    struct msh *msh;

    struct mat_cut mat;
    struct obj_cut obj;
    struct bnd_cut bnd;

    struct val_cut ext;

    struct cnd_bnd_cut cnd_bnd;
    struct cnd_ini_cut cnd_ini;
};

int sim_new(struct sim *sim, const char *sif);

#endif // NUMX_PDE_SIM_H
