#ifndef NUMX_PDE_SIM_H
#define NUMX_PDE_SIM_H

#include <numx/pde/cnd.h>
#include <numx/pde/fem.h>
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
        enum
        {
            SIM_EXP_GNS,
        } mod;

        char dir[128];
        char pfx[64];

        struct icap ini;
        struct icap put;
    } exp;

    struct
    {
        int num;
        int hop;
    } tdd;
};

struct sim
{
    enum
    {
        SIM_ELL,
        SIM_PBC,
        SIM_HYP,
    } mod;

    struct sim_ops ops;

    struct msh *msh;
    struct fem *fem;

    struct mat_cut mat;
    struct val_cut src;
    struct obj_cut obj;
    struct bnd_cut bnd;

    struct cnd_ini_cut cnd_ini;
    struct cnd_bnd_cut cnd_bnd;
};

int sim_new(struct sim *sim);
int sim_cls(struct sim *sim);

int sim_imp_gns(struct sim *sim, const char *gns);
int sim_imp_elm(struct sim *sim, const char *sif);

int sim_exp_ini_vtu(void *ctx, int n, ...);
int sim_exp_ini_gns(void *ctx, int n, ...);

int sim_exp_put_vtu(void *ctx, int n, ...);
int sim_exp_put_gns(void *ctx, int n, ...);

int sim_run(struct sim *sim);

#endif // NUMX_PDE_SIM_H
