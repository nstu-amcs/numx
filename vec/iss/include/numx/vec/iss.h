#ifndef NUMX_VEC_ISS_H
#define NUMX_VEC_ISS_H

#include <numx/com/cap.h>
#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>

typedef enum iss_mod
{
    ISS_JAC,
    ISS_RLX,
    ISS_GMR,
    ISS_BCG,
} iss_mod;

typedef enum iss_con
{
    ISS_CON_NON,
    ISS_CON_ILU,
} iss_con;

typedef struct iss_ops
{
    enum iss_con con; // Preconditioner;
    bool         pet;

    int    max; // Maximum number of iterations.
    double err; // Convergence tolerance.

    /**
     *  Iteration callback (user defined).
     *  Calls on each iteration.
     */
    struct
    {
        void *ctx;
        void (*run)(void *ctx, struct iss_ops *ops);
    } itr;

    /**
     *  Runtime data made available by solver.
     *  Updates on each iteration.
     */
    struct
    {
        int    itr; // number of iterations
        double err; // residual
    } run;
} iss_ops;

typedef struct iss_jac_ops
{
    struct iss_ops ops;

    double rlx; // relaxation parameter
} iss_jac_ops;

typedef struct iss_rlx_ops
{
    struct iss_ops ops;

    double rlx; // relaxation parameter
} iss_rlx_ops;

typedef struct iss_gmr_ops
{
    struct iss_ops ops;

    double rst; // Iteration to restart.
} iss_gmr_ops;

typedef struct iss_bcg_ops
{
    struct iss_ops ops;

    /** Preconditioner (maybe null). */
    union con
    {
        struct imtx *im;
        struct smtx *sm;
        struct dmtx *dm;
    } con;
} iss_bcg_ops;

typedef struct iss_pet_ops
{
    union
    {
        struct iss_gmr_ops *gmr;
        struct iss_bcg_ops *bcg;
    } ops;
} iss_pet_ops;

int iiss_jac_slv(struct imtx *m, struct vec *x, struct vec *f, struct iss_jac_ops *ops);
int iiss_rlx_slv(struct imtx *m, struct vec *x, struct vec *f, struct iss_rlx_ops *ops);
int iiss_gmr_slv(struct imtx *m, struct vec *x, struct vec *f, struct iss_gmr_ops *ops);
int iiss_bcg_slv(struct imtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *ops);

int diss_jac_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_jac_ops *ops);
int diss_rlx_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_rlx_ops *ops);
int diss_gmr_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_gmr_ops *ops);
int diss_bcg_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *ops);

int siss_jac_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_jac_ops *ops);
int siss_rlx_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_rlx_ops *ops);
int siss_gmr_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_gmr_ops *ops);
int siss_bcg_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *ops);

/**
 *  Jacobi solver for systems of linear equations.
 */
#define iss_jac_slv(M, x, f, o)                                                                    \
    _Generic((M),                                                                                  \
        struct imtx *: iiss_jac_slv,                                                               \
        struct dmtx *: diss_jac_slv,                                                               \
        struct smtx *: siss_jac_slv)(M, x, f, o)

/**
 *  Relaxation (Gauss-Seidel) solver for systems of linear equations.
 */
#define iss_rlx_slv(M, x, f, o)                                                                    \
    _Generic((M),                                                                                  \
        struct imtx *: iiss_rlx_slv,                                                               \
        struct dmtx *: diss_rlx_slv,                                                               \
        struct smtx *: siss_rlx_slv)(M, x, f, o)

/**
 *  Generalized minimal residual solver for systems of linear equations.
 */
#define iss_gmr_slv(M, x, f, o)                                                                    \
    _Generic((M),                                                                                  \
        struct imtx *: iiss_gmr_slv,                                                               \
        struct dmtx *: diss_gmr_slv,                                                               \
        struct smtx *: siss_gmr_slv)(M, x, f, o)

/**
 *  Stabilized biconjugate gradient solver for systems of linear equations.
 */
#define iss_bcg_slv(M, x, f, o)                                                                    \
    _Generic((M),                                                                                  \
        struct imtx *: iiss_bcg_slv,                                                               \
        struct dmtx *: diss_bcg_slv,                                                               \
        struct smtx *: siss_bcg_slv)(M, x, f, o)

#endif // NUMX_VEC_ISS_H
