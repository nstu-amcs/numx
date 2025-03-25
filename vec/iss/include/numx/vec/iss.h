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

typedef struct iss_itr
{
    double err;
    int    num;
} iss_itr;

typedef struct iss_ops
{
    /**
     * Iteration callback closure.
     *
     * @param i current iteration
     * @param r current residual
     */
    struct vcap itr;

    double err; // target residual
    int    max; // maximum number of iterations
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

    union con
    {
        struct smtx *sm;
    } con; // preconditioner (maybe null)
} iss_bcg_ops;

void iss_itr_cap(void *ctx, int n, ...);

int diss_jac_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_jac_ops o);
int diss_rlx_slv(struct dmtx *m, struct vec *x, struct vec *f, struct iss_rlx_ops o);

int siss_bcg_slv(struct smtx *m, struct vec *x, struct vec *f, struct iss_bcg_ops o);

/** Jacobi solver for systems of linear equations, Mx = f */
#define iss_jac_slv(M, x, f, o) _Generic((M), struct dmtx *: diss_jac_slv)(M, x, f, o)

/** Relaxation (Gauss-Seidel) solver for systems of linear equations, Mx = f */
#define iss_rlx_slv(M, x, f, o) _Generic((M), struct dmtx *: diss_rlx_slv)(M, x, f, o)

/** Stabilized biconjugate gradient solver for systems of linear equations, Mx = f */
#define iss_bcg_slv(M, x, f, o) _Generic((M), struct smtx *: siss_bcg_slv)(M, x, f, o)

#endif // NUMX_VEC_ISS_H
