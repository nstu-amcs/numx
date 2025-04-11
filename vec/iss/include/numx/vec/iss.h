#ifndef NUMX_VEC_ISS_H
#define NUMX_VEC_ISS_H

#include <numx/com/cap.h>
#include <numx/vec/mtx.h>
#include <numx/vec/vec.h>

typedef enum iss_mod
{
    ISS_JAC,
    ISS_RLX,
    ISS_BCG,
} iss_mod;

typedef struct iss_ops
{
    /**
     * Iteration callback (user defined).
     *
     * @param ops - that options
     */
    struct
    {
        void *ctx;
        void (*run)(void *ctx, struct iss_ops *ops);
    } itr;

    double err; // target residual
    int    max; // maximum number of iterations

    /** Runtime data made available by solver. */
    struct
    {
        int    itr;
        double err;
    } rt;
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

typedef struct iss_bcg_ops
{
    struct iss_ops ops;

    /** Preconditioner (maybe null). */
    union con
    {
        struct smtx *sm;
    } con;
} iss_bcg_ops;

int diss_jac_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_jac_ops *ops);
int diss_rlx_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_rlx_ops *ops);
int siss_bcg_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops *ops);

/** Jacobi solver for systems of linear equations, Mx = f */
#define iss_jac_slv(M, x, f, o) _Generic((M), struct dmtx *: diss_jac_slv)(M, x, f, o)

/** Relaxation (Gauss-Seidel) solver for systems of linear equations, Mx = f */
#define iss_rlx_slv(M, x, f, o) _Generic((M), struct dmtx *: diss_rlx_slv)(M, x, f, o)

/** Stabilized biconjugate gradient solver for systems of linear equations, Mx = f */
#define iss_bcg_slv(M, x, f, o) _Generic((M), struct smtx *: siss_bcg_slv)(M, x, f, o)

#endif // NUMX_VEC_ISS_H
