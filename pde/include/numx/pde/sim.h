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
    int src;
} obj;

typedef struct bnd
{
    int cnd;
} bnd;

cut_def(mat_cut, mat);
cut_def(obj_cut, obj);
cut_def(bnd_cut, bnd);

struct sim_ops
{
    void *usr;

    struct
    {
        char dir[128];
        char pfx[64];
    } msh;

    struct
    {
        enum
        {
            SIM_EXP_GNS,
            SIM_EXP_VTU,
        } mod;

        char dir[128];
        char pfx[64];
    } exp;

    struct
    {
        int    max;
        double err;
        double rlx;
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

struct sim_ell_ops
{
    struct sim_ops ops;
};

struct sim_pbc_ops
{
    struct sim_ops ops;

    int num;
    int hop;
};

struct sim_hyp_ops
{
    struct sim_ops ops;

    int num;
    int hop;
};

struct sim
{
    enum
    {
        SIM_ELL,
        SIM_PBC,
        SIM_HYP,
    } mod;

    union
    {
        struct sim_ell_ops ell;
        struct sim_pbc_ops pbc;
        struct sim_hyp_ops hyp;
    } ops;

    struct msh *msh;

    struct mat_cut mat;
    struct val_cut src;
    struct obj_cut obj;
    struct bnd_cut bnd;

    struct cnd_ini_cut cnd_ini;
    struct cnd_bnd_cut cnd_bnd;
};

int sim_new(struct sim *sim);
int sim_cls(struct sim *sim);

int sim_imp_elm(struct sim *sim, const char *sif);

#endif // NUMX_PDE_SIM_H
